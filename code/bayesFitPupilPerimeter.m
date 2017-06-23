function [ellipseFitData] = bayesFitPupilPerimeter(perimeterVideoFileName, varargin)
% [foo] = bayesFitPupilPerimeter(perimeterVideoFileName, varargin)
%
% This routine fits an ellipse to each frame of a video that contains the
% perimeter of the ellipseFitData. First, an ellipse is fit to each frame of the
% video using a non-linear search routine, with some constraints on size
% and aspect ratio of the solution. An estimate of the standard deviation
% of the parameters of the best fitting ellipse is stored as well.
% Next, a smoothing operation is conducted which estimates the posterior
% parameter values, using the measured ellipse as the likelihood and a
% non-causal, exponentially weighted window of surrounding parameter values
% as a prior.
%
% A note on ellipse parameterization: an ellipse can be specified in
% multiple forms. We adopt the standard "explicit" form for saving the
% ellipse fit results to a file. Within the context of this routine,
% however, ellipses are considered in "transparent" form (our coinage):
%
%   center (cx,cy), area (a), eccentricity (e), angle of tilt (theta)
%
% We use this parameterization to allow us to constrain fits with regard to
% these values (specifically area and eccentricity).
%
% Inputs:
%   perimeterVideoFileName: full path to a .avi file. Points on the
%   boundary of the pupil should have a value of unity, and the frame
%   should be otherwise zero filled. A frame that has no information
%   regarding the pupil (e.g., during a blink) should be zero-filled.
%   pupilFileName: full path to the .mat file in which to save pupil
%   tracking information.
%
% Optional key/value pairs
%  'verbosity' - level of verbosity. [none, full]
%  'display' - controls a display of the fitting outcome. [none, full]
%  'ellipseTransparentLB/UB' - Define the hard upper and lower boundaries
%     for the ellipse fit, in units of pixels of the video. The center
%     points should be constrained to the size of the video. Eccentricity
%     is related to ratio of the semimajor and semiminor axes, and can be
%     calculated using:
%           eccentricity = axes2ecc(semimajor, semiminor)
%     If we wish to prevent ellipses with an aspect ratio greater than
%     1.2 : 1, this gives us an eccentricity UB threshold of ~0.55.
%  'exponentialTauParams' - The time constant (in video frames) of the
%     decaying exponential weighting functions that are used to construct
%     the non-causal prior arround each frame. A different time constant
%     can (and should) be specified for the different ellipse parameters.
%     This is because pupil poisition can change rapidly due to saccades,
%     but pupil area is expected to change slowly.

% Outputs:
%   foo: bar

%% Parse vargin for options passed here
%
% Setting 'KeepUmatched' to true means that we can pass the varargin{:})
% along from a calling routine without an error here, if the key/value
% pairs recognized by the calling routine are not needed here.
p = inputParser;
p.addRequired('perimeterVideoFileName',@ischar);
p.addParameter('ellipseFitDataFileName',[],@ischar);
p.addParameter('verbosity','none',@ischar);
p.addParameter('display','none',@ischar);
p.addParameter('finalFitVideoOutFileName',[],@ischar);
p.addParameter('videoOutFrameRate',30,@isnumeric);
p.addParameter('forceNumFrames',[],@isnumeric);
p.addParameter('ellipseTransparentLB',[0, 0, 1000, 0, -0.5*pi],@isnumeric);
p.addParameter('ellipseTransparentUB',[240,320,10000,0.42, 0.5*pi],@isnumeric);
p.addParameter('exponentialTauParams',[1, 1, 20, 5, 5],@isnumeric);
p.addParameter('constrainEccen_x_Theta',0.30,@isnumeric);
p.parse(perimeterVideoFileName,varargin{:});

%% Sanity check the parameters

nEllipseParams=5; % 5 params in the transparent ellipse form

if length(p.Results.ellipseTransparentLB)~=nEllipseParams
    error('Specify the ellipseTransparentLB with 5 parameters in transparent form');
end
if length(p.Results.ellipseTransparentUB)~=nEllipseParams
    error('Specify the ellipseTransparentUB with 5 parameters in transparent form');
end
if length(p.Results.exponentialTauParams)~=nEllipseParams
    error('Specify the exponentialTauParams with 5 parameters in transparent form');
end
if sum(p.Results.ellipseTransparentUB>=p.Results.ellipseTransparentLB)~=nEllipseParams
    error('Lower bounds must be equal to or less than upper bounds');
end


%% Conduct an initial ellipse fit for each video frame

% Alert the user
if strcmp(p.Results.verbosity,'full')
    disp('Conducting initial ellipse fit to pupil perimeter video');
end

% Create the video object for reading, check the number of frames
inVideoObj = VideoReader(perimeterVideoFileName);
numFrames = floor(inVideoObj.Duration*inVideoObj.FrameRate);

% Constrain the number of frames if requested
if ~isempty(p.Results.forceNumFrames)
    numFrames = p.Results.forceNumFrames;
end

% Create a figure to hold the fit result movie frames, and display if
% requested
if strcmp(p.Results.display,'full')
    frameFig = figure( 'Visible', 'on');
else
    frameFig = figure( 'Visible', 'off');
end

% Initialize the pupil struct
ellipseFitData.X = nan(numFrames,1);
ellipseFitData.Y = nan(numFrames,1);
ellipseFitData.area = nan(numFrames,1);
ellipseFitData.pInitialFitTransparent = nan(numFrames,5);
ellipseFitData.pInitialFitSD = nan(numFrames,5);
ellipseFitData.pFinalFitTransparent = nan(numFrames,5);
ellipseFitData.pFitExplicit = nan(numFrames,5);
ellipseFitData.FinalFitError = nan(numFrames,1);

% Create a non-linear constraint for the ellipse fit. If no parameters are
%  given, then this is an identity function that does not provide any
%  constraint

if isempty(p.Results.constrainEccen_x_Theta)
    nonlinconst = [];
else
    parameterfun = @(x) p.Results.constrainEccen_x_Theta;
    nonlinconst = @(x) restrictEccenByTheta(x,parameterfun());
end

% Loop through the frames
for ii = 1:numFrames
    
    % readFrame loads frames sequentially as it is calledl; make gray
    thisFrame = rgb2gray(readFrame(inVideoObj));
    
    % get the boundary points
    [Yc, Xc] = ind2sub(size(thisFrame),find(thisFrame));
    
    % fit an ellipse to the boundaary
    if isempty(Xc) || isempty(Yc)
        pInitialFitTransparent=[NaN,NaN,NaN,NaN,NaN];
        pFitSD=[NaN,NaN,NaN,NaN,NaN];
    else
        [pInitialFitTransparent, pFitSD, ~] = ...
            calcPupilLikelihood(Xc, Yc, ...
            p.Results.ellipseTransparentLB, ...
            p.Results.ellipseTransparentUB, ...
            nonlinconst);
    end
    
    % store results
    ellipseFitData.pInitialFitTransparent(ii,:) = pInitialFitTransparent';
    ellipseFitData.pInitialFitSD(ii,:) = pFitSD';
    
    % Plot the pupil boundary data points
    imshow(thisFrame)
    if ~isnan(pInitialFitTransparent(1))
        pFitImplicit = ellipse_ex2im(ellipse_transparent2ex(pInitialFitTransparent));
        a = num2str(pFitImplicit(1));
        b = num2str(pFitImplicit(2));
        c = num2str(pFitImplicit(3));
        d = num2str(pFitImplicit(4));
        e = num2str(pFitImplicit(5));
        f = num2str(pFitImplicit(6));
        
        % ellipse impicit equation
        eqt= ['(',a, ')*x^2 + (',b,')*x*y + (',c,')*y^2 + (',d,')*x+ (',e,')*y + (',f,')'];
        
        hold on
        h= ezplot(eqt,[1, 240, 1, 320]);
        set (h, 'Color', 'green')
    end % check for a valid ellipse fit to plot
    
end % loop over frames

% close the video object
clear RGB inVideoObj

% close the figure
close(frameFig)

%% Conduct a Bayesian smoothing operation

% Alert the user
if strcmp(p.Results.verbosity,'full')
    disp('Conducting Bayesian smoothing of fit values');
end

% Set up the decayng exponential weighting functions
% Define tau prior values for each parameter
window=max(p.Results.exponentialTauParams)*4;
windowSupport=1:1:window;
for jj=1:length(p.Results.exponentialTauParams)
    baseExpFunc=exp(-1/p.Results.exponentialTauParams(jj)*windowSupport);
    exponentialWeights(jj,:)=[fliplr(baseExpFunc) NaN baseExpFunc];
end

% Create the video object for reading
inVideoObj = VideoReader(perimeterVideoFileName);

% Create a figure to hold the fit result movie frames, and display if
% requested
if strcmp(p.Results.display,'full')
    frameFig = figure( 'Visible', 'on');
else
    frameFig = figure( 'Visible', 'off');
end

% Create the video object for writing, if requested
if ~isempty(p.Results.finalFitVideoOutFileName)
    outVideoObj = VideoWriter(p.Results.finalFitVideoOutFileName);
    outVideoObj.FrameRate = p.Results.videoOutFrameRate;
    open(outVideoObj);
end

for ii = 1:numFrames
    
    % readFrame loads frames sequentially as it is called; make gray
    thisFrame = rgb2gray(readFrame(inVideoObj));
    
    % get the boundary points
    [Yc, Xc] = ind2sub(size(thisFrame),find(thisFrame));
    
    % if this frame has no data, don't attempt to fit, else fit
    if isempty(Xc) || isempty(Yc)
        pFinalFitTransparent=[NaN,NaN,NaN,NaN,NaN];
        fitError=NaN;
    else
        % Calculate the prior, which is the mean of the surrounding fit
        % values, weighted by a decaying exponential in time
        rangeLowSignal=max([ii-window,1]);
        rangeHiSignal=min([ii+window,numFrames]);
        restrictLowWindow= max([(ii-window-1)*-1,0]);
        restrictHiWindow = max([(numFrames-ii-window)*-1,0]);
        for jj=1:5
            dataVector=squeeze(ellipseFitData.pInitialFitTransparent(:,jj))';
            dataVector=dataVector(rangeLowSignal:rangeHiSignal);
            weightFunction=exponentialWeights(jj,1+restrictLowWindow:end-restrictHiWindow);
            pPriorMeanTransparent(jj) = nansum(weightFunction.*dataVector,2)./nansum(weightFunction,2);
            pPriorSDTransparent(jj) = nanstd(dataVector,weightFunction);
        end
        
        % Retrieve the initialFit for this frame
        pInitialFitTransparent = ellipseFitData.pInitialFitTransparent(ii,:);
        pFitSD = ellipseFitData.pInitialFitSD(ii,:);
        
        % calculate the posterior values for the pupil fits, given the current
        % measurement and the priors
        pPosteriorTransparent = pPriorSDTransparent.^2.*pInitialFitTransparent./(pPriorSDTransparent.^2+pFitSD.^2) + ...
            pFitSD.^2.*pPriorMeanTransparent./(pPriorSDTransparent.^2+pFitSD.^2);
        
        % re-calculate the fit, fixing the pupil size from the posterior
        lb_pinArea = p.Results.ellipseTransparentLB; lb_pinArea(3) = pPosteriorTransparent(3);
        ub_pinArea = p.Results.ellipseTransparentUB; ub_pinArea(3) = pPosteriorTransparent(3);
        [pFinalFitTransparent, ~, fitError] = calcPupilLikelihood(Xc,Yc, lb_pinArea, ub_pinArea, nonlinconst);
    end % check if there are any perimeter points to fit
    
    % store results
    ellipseFitData.pFinalFitTransparent(ii,:) = pFinalFitTransparent';
    ellipseFitData.FinalFitError(ii) = fitError;
    pFitExplicit = (ellipse_transparent2ex(pFinalFitTransparent))';
    ellipseFitData.pFitExplicit(ii,:)= pFitExplicit;
    ellipseFitData.X(ii) = pFitExplicit(1);
    ellipseFitData.Y(ii) = pFitExplicit(2);
    ellipseFitData.area(ii) = pFinalFitTransparent(3);
    
    % Plot the pupil boundary data points
    imshow(thisFrame)
    
    if ~isnan(pFinalFitTransparent(1))
        pFitImplicit = ellipse_ex2im(ellipse_transparent2ex(pFinalFitTransparent));
        a = num2str(pFitImplicit(1));
        b = num2str(pFitImplicit(2));
        c = num2str(pFitImplicit(3));
        d = num2str(pFitImplicit(4));
        e = num2str(pFitImplicit(5));
        f = num2str(pFitImplicit(6));
        
        % ellipse impicit equation
        eqt= ['(',a, ')*x^2 + (',b,')*x*y + (',c,')*y^2 + (',d,')*x+ (',e,')*y + (',f,')'];
        
        hold on
        h= ezplot(eqt,[1, 240, 1, 320]);
        set (h, 'Color', 'green')
    end
    
    % save frame
    if ~isempty(p.Results.finalFitVideoOutFileName)
        frame   = getframe(frameFig);
        writeVideo(outVideoObj,frame);
    end % check if we are saving an movie out
    
end % loop over frames to calculate the posterior

% close the figure
close(frameFig)

% close the inVideoObj
clear RGB inVideoObj

% close the outVideoObj
if ~isempty(p.Results.finalFitVideoOutFileName)
    close(outVideoObj);
end

% save the ellipse fit results if requested
if ~isempty(p.Results.ellipseFitDataFileName)
    save(p.Results.ellipseFitDataFileName,'ellipseFitData')
end

end % function


% This function implements the non-linear constraint upon the ellipse fit
% to the pupil boundary. The goal of the limit is to constrain theta to the
% cardinal axes, and more severely constrain eccentricity in the horizontal
% as compared to the vertical direction.

function [c, ceq]=restrictEccenByTheta(transparentEllipseParams,constrainEccen_x_Theta)

cardinalTolerance = (2*pi)/180;

if isempty(constrainEccen_x_Theta)
    c=[];
    ceq=[];
else
    % We implement two constraints:
    %  - the theta is on a cardinal axis (i.e., theta is from the set [-pi/2,0,pi/2])
    %  - when theta is horizontal (=0), we require that eccen be less than the
    %  more stringent horizontal eccentricity value
    c=[];
    ceq = mod(transparentEllipseParams(5),(pi/2));
    if transparentEllipseParams(5) == 0
        ceq = double(transparentEllipseParams(4) > constrainEccen_x_Theta);
    end
end

end

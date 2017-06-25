function [ellipseFitData] = bayesFitPupilPerimeter(perimeterVideoFileName, varargin)
% [foo] = bayesFitPupilPerimeter(perimeterVideoFileName, varargin)
%
% This routine fits an ellipse to each frame of a video that contains the
% perimeter of a pupil. First, an ellipse is fit to each frame of the
% video using a non-linear search routine, with some constraints on size
% and aspect ratio of the solution. An estimate of the standard deviation
% of the parameters of the best fitting ellipse is stored as well.
% Next, a smoothing operation is conducted that estimates the posterior
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
%   perimeterVideoFileName: full path to an .avi file. Points on the
%   boundary of the pupil should have a value of unity, and the frame
%   should be otherwise zero filled. A frame that has no information
%   regarding the pupil (e.g., during a blink) should be zero-filled.
%   pupilFileName: full path to the .mat file in which to save pupil
%   tracking information.
%
% Optional key/value pairs (display and I/O)
%  'verbosity' - level of verbosity. [none, full]
%  'display' - controls a display of the fitting outcome. [none, full]
%  'finalFitVideoOutFileName' - File name to save video showing the fits.
%     Defaults to empty, for which no file is saved.
%  'videoOutFrameRate' - frame rate (in Hz) of saved video. Default 30.
%  'forceNumFrames' - analyze fewer than the total number of video frames.
%
% Optional key/value pairs (analysis parameters)
%  'ellipseTransparentLB/UB' - Define the hard upper and lower boundaries
%     for the ellipse fit, in units of pixels of the video. The center
%     points should be constrained to the size of the video. Eccentricity
%     is related to ratio of the semimajor and semiminor axes, and can be
%     calculated using:
%           eccentricity = axes2ecc(semimajor, semiminor)
%     If we wish to prevent ellipses with an aspect ratio greater than
%     1.15 : 1, this gives us an eccentricity UB threshold of ~0.5.
%  'exponentialTauParams' - The time constant (in video frames) of the
%     decaying exponential weighting functions that are used to construct
%     the non-causal prior arround each frame. A different time constant
%     can (and should) be specified for the different ellipse parameters.
%     This is because pupil poisition can change rapidly due to saccades,
%     but pupil area is expected to change slowly.
%   'constrainEccen_x_Theta' - If defined, the ellipse fitting will be
%     constrained to allow only eccentric ellipses aligned with the
%     vertical or horizontal axes. Further, the two values provided will
%     limit the eccentricity of ellipses on the horizontal and vertical
%     axes, respectively.
%   'likelihoodErrorExponent' - The SD of the parameters estimated for each
%     frame are raises to this exponent, to either to weaken (>1) or
%     strengthen (<1) the influence of the current measure on the
%     posterior.
%   'nSplits' - The number of tests upon the spatial split-halves of the
%     pupil boundary values to examine to estimate a likelihood SD.
%   'nBoots' - The number of bootstrap resamples of the pupil boundary
%     points to perform to estimate a likelihood SD. This was found in
%     testing to not be useful, as the pupil boundary is oversampled, so
%     this is left at a default of zero.
%   'useParallel' - If set to true, use the Matlab parallel pool for the
%     bootstrap estimate of SD.
%   'debugMode' - If set to true, the routine attempts to load a
%     pre-existing set of initial ellipse measures (and SDs upon those
%     params), rather than re-computing these.

% Outputs:
%   ellipseFitData: A structure with multiple fields corresponding to the
%   parameters, SDs, and errors of the initial and final ellipse fits.

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
p.addParameter('ellipseTransparentUB',[240,320,10000,0.417, 0.5*pi],@isnumeric);
p.addParameter('exponentialTauParams',[0.25, 0.25, 5, 1, 1],@isnumeric);
p.addParameter('constrainEccen_x_Theta',[0.305,0.417],@isnumeric);
p.addParameter('likelihoodErrorExponent',1.25,@isnumeric);
p.addParameter('nSplits',8,@isnumeric);
p.addParameter('nBoots',0,@isnumeric);
p.addParameter('useParallel',false,@islogical);
p.addParameter('debugMode',true,@islogical);
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

%% Prepare anonymous functions that are used throughout
% Create a non-linear constraint for the ellipse fit. If no parameters are
%  given, then this is an identity function that does not provide any
%  constraint
if isempty(p.Results.constrainEccen_x_Theta)
    nonlinconst = [];
else
    nonlinconst = @(x) restrictEccenByTheta(x,p.Results.constrainEccen_x_Theta);
end

% Create an anonymous function for ellipse fitting
obtainPupilLikelihood = @(x,y) constrainedEllipseFit(x, y, ...
    p.Results.ellipseTransparentLB, p.Results.ellipseTransparentUB, nonlinconst);

% Create an anonymous function to return a rotation matrix given theta in
% radians
returnRotMat = @(theta) [cos(theta) -sin(theta); sin(theta) cos(theta)];

%% Open the video in and figure out how many frames we  will be using
% Create the video object for reading, check the number of frames
inVideoObj = VideoReader(perimeterVideoFileName);
numFrames = floor(inVideoObj.Duration*inVideoObj.FrameRate);

% Constrain the number of frames if requested
if ~isempty(p.Results.forceNumFrames)
    numFrames = p.Results.forceNumFrames;
end

%% Conduct (or load) an initial ellipse fit for each video frame
if p.Results.debugMode
    load(p.Results.ellipseFitDataFileName);
else
    % Alert the user
    if strcmp(p.Results.verbosity,'full')
        disp('Conducting initial ellipse fit to pupil perimeter video');
    end
    
    % Create a figure to hold the fit result movie frames, and display if
    % requested
    if strcmp(p.Results.display,'full')
        frameFig = figure( 'Visible', 'on');
    else
        frameFig = figure( 'Visible', 'off');
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
            pInitialFitHessianSD=[NaN,NaN,NaN,NaN,NaN];
            pInitialFitSplitsSD=[NaN,NaN,NaN,NaN,NaN];
            pInitialFitBootsSD=[NaN,NaN,NaN,NaN,NaN];
        else
            % Obtain the fit to the veridical data
            [pInitialFitTransparent, pInitialFitHessianSD, ~] = ...
                obtainPupilLikelihood(Xc, Yc);
            
            if p.Results.nSplits == 0
                pInitialFitSplitsSD=[NaN,NaN,NaN,NaN,NaN];
            else
                % Find the center of the pupil boundary points, place the boundary
                % points in a matrix and shift them to the center position
                xCenter=mean(Xc); yCenter=mean(Yc);
                centerMatrix = repmat([xCenter'; yCenter'], 1, length(Xc));
                
                % Loop through nSplits and take half the data for a random rotation
                for ss=1:p.Results.nSplits
                    theta=((pi/2)/p.Results.nSplits)*ss;
                    forwardPoints = returnRotMat(theta) * ([Xc,Yc]' - centerMatrix) + centerMatrix;
                    splitIdx1 = find((forwardPoints(1,:) < median(forwardPoints(1,:))))';
                    splitIdx2 = find((forwardPoints(1,:) >= median(forwardPoints(1,:))))';
                    
                    pFitTransparentSplit(1,ss,:) = ...
                        obtainPupilLikelihood(Xc(splitIdx1), Yc(splitIdx1));
                    pFitTransparentSplit(2,ss,:) = ...
                        obtainPupilLikelihood(Xc(splitIdx2), Yc(splitIdx2));
                end % loop through splits
                
                % Calculate the SD of the parameters across splits, scaling by
                % sqrt(2) to roughly account for our use of just half the data
                pInitialFitSplitsSD=std(reshape(pFitTransparentSplit,ss*2,nEllipseParams))/sqrt(2);
            end % check if we want to do splits
            
            % Obtain the SEM of the parameters through a bootstrap resample
            if p.Results.nBoots == 0
                pInitialFitBootsSD=[NaN,NaN,NaN,NaN,NaN];
            else
                bootOptSet = statset('UseParallel',p.Results.useParallel);
                pInitialFitBootsSD = std(bootstrp(p.Results.nBoots,obtainPupilLikelihood,Xc,Yc,'Options',bootOptSet));
            end % check if we want to do bootstraps
            
        end % check if there are pupil boundary data to be fit
        
        % store results
        ellipseFitData.pInitialFitTransparent(ii,:) = pInitialFitTransparent';
        ellipseFitData.pInitialFitHessianSD(ii,:) = pInitialFitHessianSD';
        ellipseFitData.pInitialFitSplitsSD(ii,:) = pInitialFitSplitsSD';
        ellipseFitData.pInitialFitBootsSD(ii,:) = pInitialFitBootsSD';
        
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
    
    % close the figure
    close(frameFig)
    
end % debug check

% close the video object
clear RGB inVideoObj

%% Conduct a Bayesian smoothing operation
% Alert the user
if strcmp(p.Results.verbosity,'full')
    disp('Conducting Bayesian smoothing of fit values');
end

% Set up the decaying exponential weighting functions
% Define tau prior values for each parameter
window=max(p.Results.exponentialTauParams)*8;
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
        pPriorMeanTransparent=[NaN,NaN,NaN,NaN,NaN];
        pPriorSDTransparent=[NaN,NaN,NaN,NaN,NaN];
        fitError=NaN;
    else
        % Calculate the prior. The prior mean is given by the surrounding
        % fit values, weighted by a decaying exponential in time and the
        % inverse of the standard deviation of each measure. The prior
        % standard deviation is weighted only by time.
        rangeLowSignal=max([ii-window,1]);
        rangeHiSignal=min([ii+window,numFrames]);
        restrictLowWindow= max([(ii-window-1)*-1,0]);
        restrictHiWindow = max([(numFrames-ii-window)*-1,0]);
        for jj=1:5
            dataVector=squeeze(ellipseFitData.pInitialFitTransparent(:,jj))';
            dataVector=dataVector(rangeLowSignal:rangeHiSignal);
            precisionVector=squeeze(ellipseFitData.pInitialFitSplitsSD(:,jj))';
            precisionVector=precisionVector.^(-1);
            precisionVector=precisionVector(rangeLowSignal:rangeHiSignal);
            precisionVector=precisionVector-nanmin(precisionVector);
            precisionVector=precisionVector/nanmax(precisionVector);
            temporalWeightVector=exponentialWeights(jj,1+restrictLowWindow:end-restrictHiWindow);
            combinedWeightVector=precisionVector.*temporalWeightVector;
            pPriorMeanTransparent(jj) = nansum(dataVector.*combinedWeightVector,2)./ ...
                nansum(combinedWeightVector(~isnan(dataVector)),2);
            pPriorSDTransparent(jj) = nanstd(dataVector,temporalWeightVector);
        end
        
        % Retrieve the initialFit for this frame
        pInitialFitTransparent = ellipseFitData.pInitialFitTransparent(ii,:);
        pInitialFitSplitsSD = ellipseFitData.pInitialFitSplitsSD(ii,:);
        
        % Raise the estimate of the SD from the initial fit to a
        % passed exponent. This is used to adjust the relative weighting of
        % the current frame realtive to the prior
        pInitialFitSplitsSD = pInitialFitSplitsSD .^ p.Results.likelihoodErrorExponent;
        
        % calculate the posterior values for the pupil fits, given the current
        % measurement and the priors
        pPosteriorTransparent = pPriorSDTransparent.^2.*pInitialFitTransparent./(pPriorSDTransparent.^2+pInitialFitSplitsSD.^2) + ...
            pInitialFitSplitsSD.^2.*pPriorMeanTransparent./(pPriorSDTransparent.^2+pInitialFitSplitsSD.^2);
        
        % refit the data points to deal with any nans in the posterior, and
        % to obtain a measure of the fit error
        lb_pin = p.Results.ellipseTransparentLB;
        ub_pin = p.Results.ellipseTransparentUB;
        lb_pin(~isnan(pPosteriorTransparent))=pPosteriorTransparent(~isnan(pPosteriorTransparent));
        ub_pin(~isnan(pPosteriorTransparent))=pPosteriorTransparent(~isnan(pPosteriorTransparent));
        [pFinalFitTransparent, ~, fitError] = constrainedEllipseFit(Xc,Yc, lb_pin, ub_pin, nonlinconst);

    end % check if there are any perimeter points to fit
    
    % store results
    ellipseFitData.pFinalFitTransparent(ii,:) = pFinalFitTransparent';
    ellipseFitData.finalFitError(ii) = fitError;
    ellipseFitData.pPriorMeanTransparent(ii,:)= pPriorMeanTransparent';
    ellipseFitData.pPriorSDTransparent(ii,:)= pPriorSDTransparent';
    
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
    end % check if we are saving a movie out
    
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




function [c, ceq]=restrictEccenByTheta(transparentEllipseParams,constrainEccen_x_Theta)
% function [c, ceq]=restrictEccenByTheta(transparentEllipseParams,constrainEccen_x_Theta)
%
% This function implements a non-linear constraint upon the ellipse fit
% to the pupil boundary. The goal of the limit is to constrain theta to the
% cardinal axes, and more severely constrain eccentricity in the horizontal
% as compared to the vertical direction.


% First constraint (equality)
%  - the theta is on a cardinal axis (i.e., theta is from the set [-pi/2,0,pi/2])
ceq = mod(abs(transparentEllipseParams(5)),(pi/2));

% Second constraint (inequality)
%  - require the eccen to be less than constrainEccen_x_Theta, where this
%    has one value for horizontal ellipses (i.e., abs(theta)<pi/2) and a
%    second value for vertical ellipses.
if abs(transparentEllipseParams(5))<(pi/4)
    c = transparentEllipseParams(4) - constrainEccen_x_Theta(1);
else
    c = transparentEllipseParams(4) - constrainEccen_x_Theta(2);
end

end

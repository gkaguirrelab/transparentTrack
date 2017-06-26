function [ellipseFitData] = bayesFitPupilPerimeter(perimeterVideoFileName, varargin)
% [ellipseFitData] = bayesFitPupilPerimeter(perimeterVideoFileName, varargin)
%
% This routine fits an ellipse to each frame of a video that contains the
% perimeter of the pupil. First, an ellipse is fit to each frame of the
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
% NOTES REGARDING USE OF PARALLEL POOL
%
% The initial ellipse fitting is conducted within a parfor loop. The
% parallel pool will not be used unless the key/value pair
% 'useParallel' is set to true. The routine should gracefully fall-back on
% serial processing if the parallel pool is unavailable.
%
% To use the parallel pool with TbTb, provide the identity of the repo
% name in the 'tbtbRepoName', which is then used to configure the workers.

%
% INPUTS:
%   perimeterVideoFileName: full path to an .avi file. Points on the
%     boundary of the pupil should have a value of unity, and the frame
%     should be otherwise zero-filled. A frame that has no information
%     regarding the pupil (e.g., during a blink) should be zero-filled.
%
% Optional key/value pairs (display and I/O)
%  'verbosity' - level of verbosity. [none, full]
%  'display' - controls a display of the fitting outcome. [none, full]
%  'ellipseFitDataFileName': full path to the .mat file in which to save
%     pupil tracking information.
%  'finalFitVideoOutFileName' - File name to save video showing the fits.
%     Defaults to empty, in which case no file is saved.
%  'videoOutFrameRate' - frame rate (in Hz) of saved video. Default 30.
%
% Optional key/value pairs (flow control)
% 
%  'forceNumFrames' - analyze fewer than the total number of video frames.
%  'useParallel' - If set to true, use the Matlab parallel pool for the
%    initial ellipse fitting.
%  'tbtbProjectName' - The workers in the parallel pool are configured by
%    issuing a tbUseProject command for the project specified here.
%  'developmentMode' - If set to true, the routine attempts to load a
%    pre-existing set of initial ellipse measures (and SDs upon those
%    params), rather than re-computing these. This allows more rapid
%    exploration of parameter settigns that guide the Bayesian smoothing.
%
% Optional key/value pairs (fitting parameters)
%
%  'ellipseTransparentLB/UB' - Define the hard upper and lower boundaries
%     for the ellipse fit, in units of pixels of the video. The center
%     points should be constrained to the size of the video. Eccentricity
%     is related to ratio of the semimajor and semiminor axes, and can be
%     calculated using:
%           eccentricity = axes2ecc(semimajor, semiminor)
%     If we wish to prevent ellipses with an aspect ratio greater than
%     1.1 : 1, this gives us an eccentricity UB threshold of ~0.417.
%  'exponentialTauParams' - The time constant (in video frames) of the
%     decaying exponential weighting functions that are used to construct
%     the non-causal prior arround each frame. A different time constant
%     can (and should) be specified for the different ellipse parameters.
%     This is because pupil poisition can change rapidly due to saccades,
%     but pupil area is expected to change slowly.
%   'constrainEccen_x_Theta' - If defined, the ellipse fitting will be
%     constrained to allow only eccentric ellipses aligned with the
%     vertical or horizontal axes. Further, the two values provided will
%     differently limit the eccentricity of ellipses on the horizontal and
%     vertical axes, respectively.
%   'likelihoodErrorExponent' - The SD of the parameters estimated for each
%     frame are raised to this exponent, to either to weaken (>1) or
%     strengthen (<1) the influence of the current measure on the
%     posterior.
%   'nSplits' - The number of tests upon the spatial split-halves of the
%     pupil boundary values to examine to estimate a likelihood SD.
%   'nBoots' - The number of bootstrap resamples of the pupil boundary
%     points to perform to estimate a likelihood SD. This was found in
%     testing to not be useful, as the pupil boundary is oversampled, so
%     this is left at a default of zero.
%   'priorCenterNaN' - Controls the behavior of the weighting function for
%     the prior for the current time point. See comments below for details.
%   'whichLikelihoodSD' - The variance of the measured parameters for a
%     frame can be estimated using different methods. This setting controls
%     which of these is used to set the SD of the likelihood in the
%     calculation of the posterior. Valid values:
%       'pInitialFitHessianSD'
%       'pInitialFitSplitsSD'
%       'pInitialFitBootsSD'
%
% OUTPUTS:
%   ellipseFitData: A structure with multiple fields corresponding to the
%     parameters, SDs, and errors of the initial and final ellipse fits.

%% Parse vargin for options passed here
p = inputParser;

% Required
p.addRequired('perimeterVideoFileName',@ischar);

% Optional display and I/O params
p.addParameter('verbosity','none',@ischar);
p.addParameter('display','none',@ischar);
p.addParameter('ellipseFitDataFileName',[],@ischar);
p.addParameter('finalFitVideoOutFileName',[],@ischar);
p.addParameter('videoOutFrameRate',30,@isnumeric);

% Optional flow control params
p.addParameter('forceNumFrames',[],@isnumeric);
p.addParameter('useParallel',false,@islogical);
p.addParameter('tbtbRepoName','LiveTrackAnalysisToolbox',@ischar);
p.addParameter('developmentMode',false,@islogical);

% Optional fitting params
p.addParameter('ellipseTransparentLB',[0, 0, 1000, 0, -0.5*pi],@isnumeric);
p.addParameter('ellipseTransparentUB',[240,320,10000,0.417, 0.5*pi],@isnumeric);
p.addParameter('exponentialTauParams',[.25, .25, 5, 1, 1],@isnumeric);
p.addParameter('constrainEccen_x_Theta',[0.305,0.417],@isnumeric);
p.addParameter('likelihoodErrorExponent',1.25,@isnumeric);
p.addParameter('nSplits',8,@isnumeric);
p.addParameter('nBoots',0,@isnumeric);
p.addParameter('priorCenterNaN',true,@islogical);
p.addParameter('whichLikelihoodSD','pInitialFitSplitsSD',@ischar);

% Parse the parameters
p.parse(perimeterVideoFileName,varargin{:});

%% Sanity check the parameters
nEllipseParams=5; % 5 params in the transparent ellipse form

if length(p.Results.ellipseTransparentLB)~=nEllipseParams
    error('Wrong number of elements in ellipseTransparentLB');
end
if length(p.Results.ellipseTransparentUB)~=nEllipseParams
    error('Wrong number of elements in ellipseTransparentUB');
end
if length(p.Results.exponentialTauParams)~=nEllipseParams
    error('Wrong number of elements in exponentialTauParams');
end
if sum(p.Results.ellipseTransparentUB>=p.Results.ellipseTransparentLB)~=nEllipseParams
    error('Lower bounds must be equal to or less than upper bounds');
end

%% Set up the parallel pool
if p.Results.useParallel
    poolObj = gcp;
    if isempty(poolObj)
        nWorkers=0;
    else
        nWorkers = poolObj.NumWorkers;
        % Use TbTb to configure the workers.
        if ~isempty(p.Results.tbtbRepoName)
            if strcmp(p.Results.verbosity,'full')
                fprintf('\n');
                fprintf('Configuration messages from the workers:\n');
            end
            spmd
                tbUse(p.Results.tbtbRepoName,'reset','full','verbose',false,'online',false);
            end
            if strcmp(p.Results.verbosity,'full')                
                fprintf('CAUTION: TbTb has verbose set to false.\n');
            end
        end
        % Silence warnings regarding temporary variables in the parFor loop
        spmd
            warning('off','MATLAB:mir_warning_maybe_uninitialized_temporary');
        end
    end
else
    nWorkers=0;
end

%% Prepare some anonymous functions
% Create a non-linear constraint for the ellipse fit. If no parameters are
% given, then create an empty function handle (and thus have no non-linear
% constraint)
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

%% Open and load the video
inVideoObj = VideoReader(perimeterVideoFileName);

% Get video dimensions
videoSizeX = inVideoObj.Width;
videoSizeY = inVideoObj.Height;

% Determine how many frames to process
nFrames = floor(inVideoObj.Duration*inVideoObj.FrameRate);
if ~isempty(p.Results.forceNumFrames)
    nFrames = p.Results.forceNumFrames;
end

% Load the entire pupil perimeter video into memory (about 50 MB for five
% minutes of 60 Hz data)
if strcmp(p.Results.verbosity,'full')
    tic
    fprintf('\n');
    fprintf(['Loading pupil perimeter file. Started ' char(datetime('now')) '\n']);
end
for ii = 1:nFrames
    % readFrame loads frames sequentially as it is called; make gray
    pupilBoundaryData(:,:,ii) = rgb2gray(readFrame(inVideoObj));
end
clear RGB inVideoObj
if strcmp(p.Results.verbosity,'full')
    toc
    fprintf('\n');
end


%% Load or calculate an initial ellipse fit for each video frame
if p.Results.developmentMode
    load(p.Results.ellipseFitDataFileName);
else
    
    % Alert the user
    if strcmp(p.Results.verbosity,'full')
        tic
        fprintf(['Initial ellipse fit. Started ' char(datetime('now')) '\n']);
        fprintf('| 0                      50                   100%% |\n');
        fprintf('.\n');
    end
    
    % Create a figure to display the fit results
    if ~p.Results.useParallel && strcmp(p.Results.display,'full')
        frameFig = figure( 'Visible', 'on');
    end
    
    % Loop through the frames
    parfor (ii = 1:nFrames, nWorkers)
        
        if strcmp(p.Results.verbosity,'full')
            if mod(ii,round(nFrames/50))==0
                fprintf('\b.\n');
            end
        end
        
        % get the data frame
        thisFrame = squeeze(pupilBoundaryData(:,:,ii));
        
        % get the boundary points
        [Yc, Xc] = ind2sub(size(thisFrame),find(thisFrame));
        
        % fit an ellipse to the boundary (if any points exist)
        if isempty(Xc) || isempty(Yc)
            pInitialFitTransparent=NaN(1,nEllipseParams);
            pInitialFitHessianSD=NaN(1,nEllipseParams);
            pInitialFitSplitsSD=NaN(1,nEllipseParams);
            pInitialFitBootsSD=NaN(1,nEllipseParams);
        else
            % Obtain the fit to the veridical data
            [pInitialFitTransparent, pInitialFitHessianSD, ~] = ...
                feval(obtainPupilLikelihood,Xc, Yc);
            
            % Re-calculate fit for splits of data points, if requested
            if p.Results.nSplits == 0
                pInitialFitSplitsSD=NaN(1,nEllipseParams);
            else
                % Find the center of the pupil boundary points, place the boundary
                % points in a matrix and shift them to the center position
                xCenter=mean(Xc); yCenter=mean(Yc);
                centerMatrix = repmat([xCenter'; yCenter'], 1, length(Xc));
                
                % Rotate the data and split in half through the center
                pFitTransparentSplit=NaN(2,p.Results.nSplits,nEllipseParams);
                for ss=1:p.Results.nSplits
                    theta=((pi/2)/p.Results.nSplits)*ss;
                    forwardPoints = feval(returnRotMat,theta) * ([Xc,Yc]' - centerMatrix) + centerMatrix;
                    splitIdx1 = find((forwardPoints(1,:) < median(forwardPoints(1,:))))';
                    splitIdx2 = find((forwardPoints(1,:) >= median(forwardPoints(1,:))))';
                    
                    pFitTransparentSplit(1,ss,:) = ...
                        feval(obtainPupilLikelihood,Xc(splitIdx1), Yc(splitIdx1));
                    pFitTransparentSplit(2,ss,:) = ...
                        feval(obtainPupilLikelihood,Xc(splitIdx2), Yc(splitIdx2));
                end % loop through splits
                
                % Calculate the SD of the parameters across splits, scaling by
                % sqrt(2) to roughly account for our use of just half the data
                pInitialFitSplitsSD=std(reshape(pFitTransparentSplit,ss*2,nEllipseParams))/sqrt(2);
            end % check if we want to do splits
            
            % Obtain the SD of the parameters through a bootstrap resample
            % of data points if requested
            if p.Results.nBoots == 0
                pInitialFitBootsSD=NaN(1,nEllipseParams);
            else
                bootOptSet = statset('UseParallel',p.Results.useParallel);
                pInitialFitBootsSD = std(bootstrp(p.Results.nBoots,obtainPupilLikelihood,Xc,Yc,'Options',bootOptSet));
            end % check if we want to do bootstraps
            
        end % check if there are pupil boundary data to be fit
        
        % store results
        loopVar_pInitialFitTransparent(ii,:) = pInitialFitTransparent';
        loopVar_pInitialFitHessianSD(ii,:) = pInitialFitHessianSD';
        loopVar_pInitialFitSplitsSD(ii,:) = pInitialFitSplitsSD';
        loopVar_pInitialFitBootsSD(ii,:) = pInitialFitBootsSD';
        
        % Plot the pupil boundary data
        if ~p.Results.useParallel && strcmp(p.Results.display,'full')
            imshow(thisFrame)
            if ~isnan(pInitialFitTransparent(1))
                % build ellipse impicit equation
                pFitImplicit = ellipse_ex2im(ellipse_transparent2ex(pInitialFitTransparent));
                a = num2str(pFitImplicit(1)); b = num2str(pFitImplicit(2));
                c = num2str(pFitImplicit(3)); d = num2str(pFitImplicit(4));
                e = num2str(pFitImplicit(5)); f = num2str(pFitImplicit(6));
                eqt = ['(',a, ')*x^2 + (',b,')*x*y + (',c,')*y^2 + (',d,')*x+ (',e,')*y + (',f,')'];
                
                % superimpose the ellipse using ezplot
                hold on
                h = ezplot(eqt,[1, videoSizeY, 1, videoSizeX]);
                set (h, 'Color', 'green')
            end % check for a valid ellipse fit to plot
        end % display check
        
    end % loop over frames
    
    % close the figure
    if ~p.Results.useParallel && strcmp(p.Results.display,'full')
        close(frameFig)
    end
    
    % gather the loop vars into the ellipse structure
    ellipseFitData.pInitialFitTransparent = loopVar_pInitialFitTransparent;
    ellipseFitData.pInitialFitHessianSD = loopVar_pInitialFitHessianSD;
    ellipseFitData.pInitialFitSplitsSD = loopVar_pInitialFitSplitsSD;
    ellipseFitData.pInitialFitBootsSD = loopVar_pInitialFitBootsSD;
    
    if strcmp(p.Results.verbosity,'full')
        toc
        fprintf('\n\n');
    end
    
end % developmentMode check

% save the ellipse fit results if requested
if ~isempty(p.Results.ellipseFitDataFileName)
    save(p.Results.ellipseFitDataFileName,'ellipseFitData')
end

%% Conduct a Bayesian smoothing operation

% Set up the decaying exponential weighting functions. The relatively large
% window (8 times the biggest time constant) is used to handle the case in
% which there is a stretch of missing data, in which case the long tails of
% the exponential can provide the prior.
window=max(p.Results.exponentialTauParams)*8;
if p.Results.priorCenterNaN
    windowSupport=1:1:window;
else
    windowSupport=1:1:window+1;
end
for jj=1:nEllipseParams
    baseExpFunc=exp(-1/p.Results.exponentialTauParams(jj)*windowSupport);
    
    % The weighting function is symmetric about the current time point. A
    % parameter flag switches the treatment of the current time point,
    % either excluding it and weighting most heavily the immediately
    % adjacent frames, or giving the current time point a weight of unity.
    if p.Results.priorCenterNaN
        exponentialWeights(jj,:)=[fliplr(baseExpFunc) NaN baseExpFunc];
    else
        exponentialWeights(jj,1:window+1)=fliplr(baseExpFunc);
        exponentialWeights(jj,window+1:window*2+1)=baseExpFunc;
    end
end

% Create a figure to hold the fit result movie frames, and display if
% requested
if strcmp(p.Results.display,'full')
    frameFig = figure( 'Visible', 'on');
else
    frameFig = figure( 'Visible', 'off');
end

% Alert the user
if strcmp(p.Results.verbosity,'full')
    tic
    fprintf(['Bayesian smoothing. Started ' char(datetime('now')) '\n']);
    fprintf('| 0                      50                   100%% |\n');
    fprintf('.\n');
end

% Pre-allocate a cell array to hold the finalFit video frames
finalFitVideo=cell(nFrames,1);

parfor (ii = 1:nFrames, nWorkers)
    
    % update progress
    if strcmp(p.Results.verbosity,'full')
        if mod(ii,round(nFrames/50))==0
            fprintf('\b.\n');
        end
    end
    
    % get the data frame
    thisFrame = squeeze(pupilBoundaryData(:,:,ii));
    
    % get the boundary points
    [Yc, Xc] = ind2sub(size(thisFrame),find(thisFrame));
    
    % if this frame has no data, don't attempt to fit, else fit
    if isempty(Xc) || isempty(Yc)
        pPosteriorMeanTransparent=NaN(1,nEllipseParams);
        pPosteriorSDTransparent=NaN(1,nEllipseParams);
        pPriorMeanTransparent=NaN(1,nEllipseParams);
        pPriorSDTransparent=NaN(1,nEllipseParams);
        fitError=NaN;
    else
        % Calculate the prior. The prior mean is given by the surrounding
        % fit values, weighted by a decaying exponential in time and the
        % inverse of the standard deviation of each measure. The prior
        % standard deviation is weighted only by time.
        
        % A bit of fussing with the range here to handle the start and the
        % end of the data vector
        rangeLowSignal=max([ii-window,1]);
        rangeHiSignal=min([ii+window,nFrames]);
        restrictLowWindow= max([(ii-window-1)*-1,0]);
        restrictHiWindow = max([(nFrames-ii-window)*-1,0]);
        
        for jj=1:nEllipseParams
            % Get the dataVector, restricted to the window range
            dataVector=squeeze(ellipseFitData.pInitialFitTransparent(:,jj))';
            dataVector=dataVector(rangeLowSignal:rangeHiSignal);
            
            % Build the precisionVector as the inverse of the measurement
            % SD on each frame, scaled to range within the window from zero
            % to unity. Thus, the noisiest measurement will not influence
            % the prior.
            precisionVector=squeeze(ellipseFitData.pInitialFitSplitsSD(:,jj))';
            precisionVector=precisionVector.^(-1);
            precisionVector=precisionVector(rangeLowSignal:rangeHiSignal);
            precisionVector=precisionVector-nanmin(precisionVector);
            precisionVector=precisionVector/nanmax(precisionVector);
            
            % The temporal weight vector is simply the exponential weights,
            % restricted to the available data widow
            temporalWeightVector = ...
                exponentialWeights(jj,1+restrictLowWindow:end-restrictHiWindow);
            
            % Combine the precision and time weights, and calculate the
            % prior mean
            combinedWeightVector=precisionVector.*temporalWeightVector;
            pPriorMeanTransparent(jj) = nansum(dataVector.*combinedWeightVector,2)./ ...
                nansum(combinedWeightVector(~isnan(dataVector)),2);
            
            % Obtain the standard deviation of the prior
            pPriorSDTransparent(jj) = nanstd(dataVector,temporalWeightVector);
        end
        
        % Retrieve the initialFit for this frame
        pLikelihoodMeanTransparent = ellipseFitData.pInitialFitTransparent(ii,:);
        
        % There are different measures available for the SD of the
        % parameters of the initial fit. The parameter 'whichLikelihoodSD'
        % controls which one of these is used for the likelihood
        if ~isfield(ellipseFitData,p.Results.whichLikelihoodSD)
            error('The requested estimate of fit SD is not available in ellipseFitData');
        else
            pLikelihoodSDTransparent = ellipseFitData.(p.Results.whichLikelihoodSD)(ii,:);
        end
        
        % Raise the estimate of the SD from the initial fit to an
        % exponent. This is used to adjust the relative weighting of
        % the current frame realtive to the prior
        pLikelihoodSDTransparent = pLikelihoodSDTransparent .^ p.Results.likelihoodErrorExponent;
        
        % Calculate the posterior values for the pupil fits, given the
        % likelihood and the prior
        pPosteriorMeanTransparent = pPriorSDTransparent.^2.*pLikelihoodMeanTransparent./(pPriorSDTransparent.^2+pLikelihoodSDTransparent.^2) + ...
            pLikelihoodSDTransparent.^2.*pPriorMeanTransparent./(pPriorSDTransparent.^2+pLikelihoodSDTransparent.^2);

        pPosteriorSDTransparent = sqrt((pPriorSDTransparent.^2.*pLikelihoodSDTransparent.^2) ./ ...
            (pPriorSDTransparent.^2+pLikelihoodSDTransparent.^2));
        
        % refit the data points to deal with any NaNs in the posterior and
        % to obtain a measure of the fit error
        lb_pin = p.Results.ellipseTransparentLB;
        ub_pin = p.Results.ellipseTransparentUB;
        lb_pin(~isnan(pPosteriorMeanTransparent))=pPosteriorMeanTransparent(~isnan(pPosteriorMeanTransparent));
        ub_pin(~isnan(pPosteriorMeanTransparent))=pPosteriorMeanTransparent(~isnan(pPosteriorMeanTransparent));
        [pPosteriorMeanTransparent, ~, fitError] = constrainedEllipseFit(Xc,Yc, lb_pin, ub_pin, nonlinconst);
        
    end % check if there are any perimeter points to fit
    
    % store results
    loopVar_pPosteriorMeanTransparent(ii,:) = pPosteriorMeanTransparent';
    loopVar_pPosteriorSDTransparent(ii,:) = pPosteriorSDTransparent';
    loopVar_finalFitError(ii) = fitError;
    loopVar_pPriorMeanTransparent(ii,:)= pPriorMeanTransparent';
    loopVar_pPriorSDTransparent(ii,:)= pPriorSDTransparent';
    
    % Plot the pupil boundary data points
    imshow(thisFrame)
    
    if ~isnan(pPosteriorMeanTransparent(1))
        % build ellipse impicit equation
        pFitImplicit = ellipse_ex2im(ellipse_transparent2ex(pPosteriorMeanTransparent));
        a = num2str(pFitImplicit(1)); b = num2str(pFitImplicit(2));
        c = num2str(pFitImplicit(3)); d = num2str(pFitImplicit(4));
        e = num2str(pFitImplicit(5)); f = num2str(pFitImplicit(6));
        eqt = ['(',a, ')*x^2 + (',b,')*x*y + (',c,')*y^2 + (',d,')*x+ (',e,')*y + (',f,')'];
        
        % superimpose the ellipse using ezplot
        hold on
        h = ezplot(eqt,[1, videoSizeY, 1, videoSizeX]);
        set (h, 'Color', 'green')
    end
    
    % collect the frame into the finalFitVideo
    if ~isempty(p.Results.finalFitVideoOutFileName)
        finalFitVideo{ii}   = getframe(frameFig);
    end % check if we are saving a movie out
    
end % loop over frames to calculate the posterior

% gather the loop vars into the ellipse structure
ellipseFitData.pPosteriorMeanTransparent=loopVar_pPosteriorMeanTransparent;
ellipseFitData.pPosteriorSDTransparent=loopVar_pPosteriorSDTransparent;
ellipseFitData.fitError=loopVar_finalFitError;
ellipseFitData.pPriorMeanTransparent=loopVar_pPriorMeanTransparent;
ellipseFitData.pPriorSDTransparent=loopVar_pPriorSDTransparent;

%% Cleanup and save data
% close the figure
close(frameFig)

% save the finalFitVideo
if ~isempty(p.Results.finalFitVideoOutFileName)
    % Create the video object for writing
    outVideoObj = VideoWriter(p.Results.finalFitVideoOutFileName);
    outVideoObj.FrameRate = p.Results.videoOutFrameRate;
    open(outVideoObj);
    for ii=1:nFrames
        writeVideo(outVideoObj,finalFitVideo{ii});
    end
    close(outVideoObj);
end

% save the ellipse fit results if requested
if ~isempty(p.Results.ellipseFitDataFileName)
    save(p.Results.ellipseFitDataFileName,'ellipseFitData')
end

% report completion of Bayesian analysis
if strcmp(p.Results.verbosity,'full')
    fprintf('\n');
    toc
    fprintf('\n\n');
end

% Delete the parallel pool
if p.Results.useParallel
    poolObj = gcp;
    if ~isempty(poolObj)
        delete(poolObj);
    end
end

end % function



function [c, ceq]=restrictEccenByTheta(transparentEllipseParams, constrainEccen_x_Theta)
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

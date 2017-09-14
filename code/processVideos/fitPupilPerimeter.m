function [pupilData] = fitPupilPerimeter(perimeterFileName, pupilFileName, varargin)
% fitPupilPerimeter(perimeterVideoFileName, varargin)
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
% A note on the results coordinates: matlab's intrinsic coordinates system
% is built so that the center of each pixel on the image cooresponds to the
% integer indexed position on the pixel itself, that means that a 3x3
% pixels image in intrisic coordinates will be represented in a grid with
% xlim = [0.5 3.5] and ylim = [0.5 3.5], with the origin being the top left
% corner of the image. This is done to facilitate the handling of images in
% many of the built-in image processing functions. We are interested in
% returning the results in "world units", which in our case correspond to
% the cartesian system having the origin in the top left corner of the
% frame and, referring to the above example, xlim = [0 3] and ylim = [0 3],
% with the assumption of using square pixels. To apply this conversion, the
% last stage before saving the pupil data is to subtract 0.5 pixels from
% every X transparent coordinate and add 0.5 pixel to every Y transparent
% coordinate.
% REF for intrinsic coordinates explaination: 
% https://blogs.mathworks.com/steve/2013/08/28/introduction-to-spatial-referencing/
%
% NOTES REGARDING USE OF PARALLEL POOL
%
% The initial ellipse fitting is conducted within a parfor loop. The
% parallel pool will not be used unless the key/value pair
% 'useParallel' is set to true. The routine should gracefully fall-back on
% serial processing if the parallel pool is unavailable.
%
% Each worker requires ~8 GB of memory to operate. It is important to keep
% total RAM usage below the physical memory limit to prevent swapping and
% a dramatic slow down in processing.
%
% To use the parallel pool with TbTb, provide the identity of the repo
% name in the 'tbtbRepoName', which is then used to configure the workers.
%
% INPUTS:
%   perimeterFileName: full path to a .mat file that contains the perimeter
%     data varaible. Points on the boundary of the pupil should have a
%     value of unity, and the frame should be otherwise zero-filled. A
%     frame that has no information regarding the pupil (e.g., during a
%     blink) should be zero-filled.
%   pupilFileName: full path to the .mat file in which to save
%     pupil tracking information.
%
% Optional key/value pairs (verbosity)
%  'verbosity' - level of verbosity. [none, full]
%
% Optional key/value pairs (flow control)
%
%  'nFrames' - analyze fewer than the total number of frames.
%  'useParallel' - If set to true, use the Matlab parallel pool for the
%    initial ellipse fitting.
%  'nWorkers' - Specify the number of workers in the parallel pool. If
%    undefined the default number will be used.
%  'tbtbProjectName' - The workers in the parallel pool are configured by
%    issuing a tbUseProject command for the project specified here.
%  'skipInitialPupilFit' - If set to true, the routine attempts to load a
%    pre-existing set of initial ellipse measures (and SDs upon those
%    params), rather than re-computing these. This allows more rapid
%    exploration of parameter settigns that guide the Bayesian smoothing.
%  'skipPupilBayes' - If set to true, the routine exits after the initial
%    ellipse fitting and prior to performing Bayesian smoothing 
%
% Optional key/value pairs (Environment parameters)
%  'tbSnapshot' - This should contain the output of the tbDeploymentSnapshot
%    performed upon the result of the tbUse command. This documents the
%    state of the system at the time of analysis.
%  'timestamp' - AUTOMATIC - The current time and date
%  'username' - AUTOMATIC - The user
%  'hostname' - AUTOMATIC - The host
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
% 
% OUTPUTS:
%   ellipseFitData: A structure with multiple fields corresponding to the
%     parameters, SDs, and errors of the initial and final ellipse fits.


%% Parse vargin for options passed here
p = inputParser; p.KeepUnmatched = true;

% Required
p.addRequired('perimeterFileName',@ischar);
p.addRequired('pupilFileName',@ischar);

% Optional display and I/O params
p.addParameter('verbosity','none',@ischar);

% Optional flow control params
p.addParameter('nFrames',Inf,@isnumeric);
p.addParameter('useParallel',false,@islogical);
p.addParameter('nWorkers',[],@(x)(isempty(x) | isnumeric(x)));
p.addParameter('tbtbRepoName','LiveTrackAnalysisToolbox',@ischar);
p.addParameter('skipInitialPupilFit',false,@islogical);
p.addParameter('skipPupilBayes',false,@islogical);

% Environment parameters
p.addParameter('tbSnapshot',[],@(x)(isempty(x) | isstruct(x)));
p.addParameter('timestamp',char(datetime('now')),@ischar);
p.addParameter('hostname',char(java.lang.System.getProperty('user.name')),@ischar);
p.addParameter('username',char(java.net.InetAddress.getLocalHost.getHostName),@ischar);

% Optional fitting params
p.addParameter('ellipseTransparentLB',[0, 0, 500, 0, -0.5*pi],@isnumeric);
p.addParameter('ellipseTransparentUB',[320,240,10000,0.5, 0.5*pi],@isnumeric);
p.addParameter('exponentialTauParams',[.25, .25, 5, 1, 1],@isnumeric);
p.addParameter('constrainEccen_x_Theta',[0.5,0.5],@isnumeric);
p.addParameter('likelihoodErrorExponent',1.25,@isnumeric);
p.addParameter('nSplits',8,@isnumeric);
p.addParameter('nBoots',0,@isnumeric);
p.addParameter('priorCenterNaN',true,@islogical);
p.addParameter('whichLikelihoodSD','pInitialFitSplitsSD',@ischar);


%% Parse and check the parameters
p.parse(perimeterFileName, pupilFileName, varargin{:});

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

%% Announce we are starting
if strcmp(p.Results.verbosity,'full')
    fprintf('Performing non-causal Bayesian fitting of the pupil boundary file:\n');
    fprintf(['\t' perimeterFileName '\n\n']);
end

%% Prepare some anonymous functions and load the pupil perimeter data
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

% Load the pupil perimeter data. It will be a structure variable
% "perimeter", with the fields .data and .meta
dataLoad=load(perimeterFileName);
perimeter=dataLoad.perimeter;
clear dataLoad

% determine how many frames we will process
if p.Results.nFrames == Inf
    nFrames=size(perimeter.data,3);
else
    nFrames = p.Results.nFrames;
end


%% Set up the parallel pool
if p.Results.useParallel
    if strcmp(p.Results.verbosity,'full')
        tic
        fprintf(['Opening parallel pool. Started ' char(datetime('now')) '\n']);
    end
    if isempty(p.Results.nWorkers)
        parpool;
    else
        parpool(p.Results.nWorkers);
    end
    poolObj = gcp;
    if isempty(poolObj)
        nWorkers=0;
    else
        nWorkers = poolObj.NumWorkers;
        % Use TbTb to configure the workers.
        if ~isempty(p.Results.tbtbRepoName)
            spmd
                tbUse(p.Results.tbtbRepoName,'reset','full','verbose',false,'online',false);
            end
            if strcmp(p.Results.verbosity,'full')
                fprintf('CAUTION: Any TbTb messages from the workers will not be shown.\n');
            end
        end
    end
    if strcmp(p.Results.verbosity,'full')
        toc
        fprintf('\n');
    end
else
    nWorkers=0;
end


%% Load or calculate an initial ellipse fit for each video frame
if p.Results.skipInitialPupilFit
    load(p.Results.pupilFileName);
else
    
    % Alert the user
    if strcmp(p.Results.verbosity,'full')
        tic
        fprintf(['Initial ellipse fit. Started ' char(datetime('now')) '\n']);
        fprintf('| 0                      50                   100%% |\n');
        fprintf('.\n');
    end
    
    % Loop through the frames
    parfor (ii = 1:nFrames, nWorkers)
        
        % Update progress
        if strcmp(p.Results.verbosity,'full')
            if mod(ii,round(nFrames/50))==0
                fprintf('\b.\n');
            end
        end
        try % this is to have information on which frame caused an error
            % get the data frame
            thisFrame = squeeze(perimeter.data(:,:,ii));
            
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
                    pInitialFitSplitsSD=nanstd(reshape(pFitTransparentSplit,ss*2,nEllipseParams))/sqrt(2);
                end % check if we want to do splits
                
                % Obtain the SD of the parameters through a bootstrap resample
                % of data points if requested
                if p.Results.nBoots == 0
                    pInitialFitBootsSD=NaN(1,nEllipseParams);
                else
                    bootOptSet = statset('UseParallel',p.Results.useParallel);
                    pInitialFitBootsSD = nanstd(bootstrp(p.Results.nBoots,obtainPupilLikelihood,Xc,Yc,'Options',bootOptSet));
                end % check if we want to do bootstraps
                
            end % check if there are pupil boundary data to be fit
            
            % store results
            loopVar_pInitialFitTransparent(ii,:) = pInitialFitTransparent';
            loopVar_pInitialFitHessianSD(ii,:) = pInitialFitHessianSD';
            loopVar_pInitialFitSplitsSD(ii,:) = pInitialFitSplitsSD';
            loopVar_pInitialFitBootsSD(ii,:) = pInitialFitBootsSD';
        catch ME
            warning ('Error while processing frame: %d', ii)
            rethrow(ME)
        end % try catch
    end % loop over frames
    
    % convert to world coordinates  (if this is the last stage)
    if p.Results.skipPupilBayes
        loopVar_pInitialFitTransparent(:,1) = loopVar_pInitialFitTransparent(:,1) - 0.5;
        loopVar_pInitialFitTransparent(:,2) = loopVar_pInitialFitTransparent(:,2) + 0.5;
    end
    % gather the loop vars into the ellipse structure    
    pupilData.pInitialFitTransparent = loopVar_pInitialFitTransparent;
    pupilData.pInitialFitHessianSD = loopVar_pInitialFitHessianSD;
    pupilData.pInitialFitSplitsSD = loopVar_pInitialFitSplitsSD;
    pupilData.pInitialFitBootsSD = loopVar_pInitialFitBootsSD;
    
    if strcmp(p.Results.verbosity,'full')
        toc
        fprintf('\n');
    end
    
end % skipInitialPupilFit mode check

% save the ellipse fit results if requested
if ~isempty(p.Results.pupilFileName)
    save(p.Results.pupilFileName,'pupilData')
end
    
% If we are not skipping the Bayesian fitting, proceed
if ~p.Results.skipPupilBayes
    
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
    
    % Alert the user
    if strcmp(p.Results.verbosity,'full')
        tic
        fprintf(['Bayesian smoothing. Started ' char(datetime('now')) '\n']);
        fprintf('| 0                      50                   100%% |\n');
        fprintf('.\n');
    end
    
    parfor (ii = 1:nFrames, nWorkers)
        
        % update progress
        if strcmp(p.Results.verbosity,'full')
            if mod(ii,round(nFrames/50))==0
                fprintf('\b.\n');
            end
        end
        
        % get the data frame
        thisFrame = squeeze(perimeter.data(:,:,ii));
        
        % get the boundary points
        [Yc, Xc] = ind2sub(size(thisFrame),find(thisFrame));
        
        pPosteriorMeanTransparent=NaN(1,nEllipseParams);
        pPosteriorSDTransparent=NaN(1,nEllipseParams);
        pPriorMeanTransparent=NaN(1,nEllipseParams);
        pPriorSDTransparent=NaN(1,nEllipseParams);
        pLikelihoodMeanTransparent=NaN(1,nEllipseParams);
        pLikelihoodSDTransparent=NaN(1,nEllipseParams);
        fitError=NaN;
        
        % if this frame has data, and the initial ellipse fit is not nan,
        % then proceed to calculate the posterior
        if ~isempty(Xc) &&  ~isempty(Yc) && sum(isnan(pupilData.pInitialFitTransparent(ii,:)))==0
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
                dataVector=squeeze(pupilData.pInitialFitTransparent(:,jj))';
                dataVector=dataVector(rangeLowSignal:rangeHiSignal);
                
                % Build the precisionVector as the inverse of the measurement
                % SD on each frame, scaled to range within the window from zero
                % to unity. Thus, the noisiest measurement will not influence
                % the prior.
                precisionVector=squeeze(pupilData.pInitialFitSplitsSD(:,jj))';
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
            pLikelihoodMeanTransparent = pupilData.pInitialFitTransparent(ii,:);
            
            % There are different measures available for the SD of the
            % parameters of the initial fit. The parameter 'whichLikelihoodSD'
            % controls which one of these is used for the likelihood
            if ~isfield(pupilData,p.Results.whichLikelihoodSD)
                error('The requested estimate of fit SD is not available in ellipseFitData');
            else
                pLikelihoodSDTransparent = pupilData.(p.Results.whichLikelihoodSD)(ii,:);
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
        
    end % loop over frames to calculate the posterior
    
    
    %% convert to world coordinates
        % convert to world coordinates
        loopVar_pPriorMeanTransparent(:,1) = loopVar_pPriorMeanTransparent(:,1) - 0.5;
        loopVar_pPriorMeanTransparent(:,2) = loopVar_pPriorMeanTransparent(:,2) + 0.5;
        loopVar_pPosteriorMeanTransparent(:,1) = loopVar_pPosteriorMeanTransparent(:,1) - 0.5;
        loopVar_pPosteriorMeanTransparent(:,2) = loopVar_pPosteriorMeanTransparent(:,2) + 0.5;
        if ~ p.Results.skipInitialPupilFit
            pupilData.pInitialFitTransparent(:,1) = pupilData.pInitialFitTransparent(:,1) - 0.5;
            pupilData.pInitialFitTransparent(:,2) = pupilData.pInitialFitTransparent(:,2) + 0.5;
        end
    
    %% Clean up and save the fit results
    
    % gather the loop vars into the ellipse structure
    pupilData.pPriorMeanTransparent=loopVar_pPriorMeanTransparent;
    pupilData.pPriorSDTransparent=loopVar_pPriorSDTransparent;
    pupilData.pPosteriorMeanTransparent=loopVar_pPosteriorMeanTransparent;
    pupilData.pPosteriorSDTransparent=loopVar_pPosteriorSDTransparent;
    pupilData.fitError=loopVar_finalFitError';
    
    % add a meta field with analysis details
    pupilData.meta = p.Results;
    pupilData.meta.coordinatesSystem = 'worldCoordinates';
    
    % save the ellipse fit results if requested
    if ~isempty(p.Results.pupilFileName)
        save(p.Results.pupilFileName,'pupilData')
    end
    
    % report completion of Bayesian analysis
    if strcmp(p.Results.verbosity,'full')
        toc
        fprintf('\n');
    end
    
end % check if we are skipping Bayesian smoothing


%% Delete the parallel pool
if strcmp(p.Results.verbosity,'full')
    tic
    fprintf(['Closing parallel pool. Started ' char(datetime('now')) '\n']);
end
if p.Results.useParallel
    poolObj = gcp;
    if ~isempty(poolObj)
        delete(poolObj);
    end
end
if strcmp(p.Results.verbosity,'full')
    toc
    fprintf('\n');
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

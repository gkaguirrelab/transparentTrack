function [pupilData] = smoothPupilParameters(perimeterFileName, pupilFileName, varargin)
% smoothPupilParameters(pupilFileName, varargin)
%
% This routine implements a smoothing operation that estimates posterior parameter
% values, using the measured ellipse as the likelihood and a non-causal,
% exponentially weighted window of surrounding parameter values as a prior.
% This is an empirical Bayes approach.
%
% Optionally, additional iterations of empirical Bayes smoothing can be
% conducted. Default settings cause the temporal domain of the smoothing to
% become shorter and shorter, preventing the iterative process from
% converging on a flat line.
%
% A note on ellipse parameterization: an ellipse can be specified in
% multiple forms. Within the context of this routine, and in saved files,
% ellipses are considered in "transparent" form (our coinage):
%
%   center (cx,cy), area (a), eccentricity (e), angle of tilt (theta)
%
% We use this parameterization to allow us to constrain fits with regard to
% these values (specifically area and eccentricity).
%
%
% NOTES REGARDING USE OF PARALLEL POOL
%
% The parallel pool will not be used unless the key/value pair
% 'useParallel' is set to true. The routine should gracefully fall-back on
% serial processing if the parallel pool is unavailable.
%
% Each worker requires ~8 GB of memory to operate. It is important to keep
% total RAM usage below the physical memory limit to prevent swapping and a
% dramatic slow down in processing.
%
% To use the parallel pool with TbTb, provide the identity of the repo name
% in the 'tbtbRepoName', which is then used to configure the workers.
%
% INPUTS:
%   perimeterFileName: full path to a .mat file that contains the perimeter
%     data varaible. Points on the boundary of the pupil should have a
%     value of unity, and the frame should be otherwise zero-filled. A
%     frame that has no information regarding the pupil (e.g., during a
%     blink) should be zero-filled.
%   pupilFileName: full path to the .mat file that contains the pupil data
%    to be smoothed. This file will be over-written by the output.
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
%   'nAdditionalBayes' - The number of additional passes of empirical Bayes
%     that is performed upon the posterior.
%   'shrinkTauParamFactor' - The factor by which the exponentialTauParams
%     are reduced on each iteration of the additional Bayes smoothing.
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
%   pupilData: A structure with multiple fields corresponding to the
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
p.addParameter('tbtbRepoName','transparentTrack',@ischar);

% Environment parameters
p.addParameter('tbSnapshot',[],@(x)(isempty(x) | isstruct(x)));
p.addParameter('timestamp',char(datetime('now')),@ischar);
p.addParameter('hostname',char(java.lang.System.getProperty('user.name')),@ischar);
p.addParameter('username',char(java.net.InetAddress.getLocalHost.getHostName),@ischar);

% Optional fitting params
p.addParameter('ellipseTransparentLB',[0, 0, 400, 0, -0.5*pi],@isnumeric);
p.addParameter('ellipseTransparentUB',[320,240,10000,0.5, 0.5*pi],@isnumeric);
p.addParameter('exponentialTauParams',[.25, .25, 5, 1, 1],@isnumeric);
p.addParameter('constrainEccen_x_Theta',[0.5,0.5],@isnumeric);
p.addParameter('likelihoodErrorExponent',1.25,@isnumeric);
p.addParameter('nAdditionalBayes',0,@isnumeric);
p.addParameter('shrinkTauParamFactor',4,@isnumeric);
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

% Load the pupil perimeter data. It will be a structure variable
% "perimeter", with the fields .data and .meta
dataLoad=load(perimeterFileName);
perimeter=dataLoad.perimeter;
clear dataLoad

% Load the pupil data. It will be a structure variable "pupilData"
dataLoad=load(pupilFileName);
pupilData=dataLoad.pupilData;
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


%% Conduct empirical Bayes smoothing

% There are different measures available for the SD of the
% parameters of the initial fit. The parameter 'whichLikelihoodSD'
% controls which one of these is used for the likelihood. check that
% the requested SD measure is available
if ~isfield(pupilData,p.Results.whichLikelihoodSD)
    error('The requested estimate of fit SD is not available in ellipseFitData');
end

% Set up the decaying exponential weighting functions. The relatively large
% window (8 times the biggest time constant) is used to handle the case in
% which there is a stretch of missing data, in which case the long tails of
% the exponential can provide the prior.
window=ceil(max([max(p.Results.exponentialTauParams)*8,8]));
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
    
    % initialize some variables
    pPosteriorMeanTransparent=NaN(1,nEllipseParams);
    pPosteriorSDTransparent=NaN(1,nEllipseParams);
    pPriorMeanTransparent=NaN(1,nEllipseParams);
    pPriorSDTransparent=NaN(1,nEllipseParams);
    pLikelihoodMeanTransparent=NaN(1,nEllipseParams);
    pLikelihoodSDTransparent=NaN(1,nEllipseParams);
    fitError=NaN;
    
    % get the data frame
    thisFrame = squeeze(perimeter.data(:,:,ii));
    
    % get the boundary points
    [Yc, Xc] = ind2sub(size(thisFrame),find(thisFrame));
    
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
            precisionVector=squeeze(pupilData.(p.Results.whichLikelihoodSD)(:,jj))';
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
        pLikelihoodSDTransparent = pupilData.(p.Results.whichLikelihoodSD)(ii,:);
        
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
        [pPosteriorMeanTransparent, ~, fitError] = constrainedEllipseFit(Xc,Yc, lb_pin, ub_pin, []);
        
    end % check if there are any perimeter points to fit
    
    % store results
    loopVar_pPosteriorMeanTransparent(ii,:) = pPosteriorMeanTransparent';
    loopVar_pPosteriorSDTransparent(ii,:) = pPosteriorSDTransparent';
    loopVar_finalFitError(ii) = fitError;
    loopVar_pPriorMeanTransparent(ii,:)= pPriorMeanTransparent';
    loopVar_pPriorSDTransparent(ii,:)= pPriorSDTransparent';
    
end % loop over frames to calculate the posterior

% Conduct additional passes of Bayes smoothing if requested
if p.Results.nAdditionalBayes > 0
    for pp = 1:p.Results.nAdditionalBayes
        
        % Copy the posterior from the last loop into the likelihood for
        % this loop
        loopVar_pLikelihoodMeanFromLastLoop = loopVar_pPosteriorMeanTransparent;
        loopVar_pLikelihoodSDFromLastLoop = loopVar_pPosteriorSDTransparent;
        
        % Build a new set of exponential decay weights, with the tau
        % parameters reduced as the loops progress
        clear exponentialWeights
        thisLoopTauParams = p.Results.exponentialTauParams ./ (p.Results.shrinkTauParamFactor.^pp);
        window=ceil(max([max(thisLoopTauParams)*8,8]));
        if p.Results.priorCenterNaN
            windowSupport=1:1:window;
        else
            windowSupport=1:1:window+1;
        end
        for jj=1:nEllipseParams
            baseExpFunc=exp(-1/thisLoopTauParams(jj)*windowSupport);
            
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
        
        % update progress
        if strcmp(p.Results.verbosity,'full')
            fprintf(['Additional Bayes pass ' num2str(pp) '\n\n']);
        end
        
        parfor (ii = 1:nFrames, nWorkers)
            
            % update progress
            if strcmp(p.Results.verbosity,'full')
                if mod(ii,round(nFrames/50))==0
                    fprintf('\b.\n');
                end
            end % check for verbosity
            
            % initialize some variables
            pPosteriorMeanTransparent=NaN(1,nEllipseParams);
            pPosteriorSDTransparent=NaN(1,nEllipseParams);
            pPriorMeanTransparent=NaN(1,nEllipseParams);
            pPriorSDTransparent=NaN(1,nEllipseParams);
            fitError=NaN;
            
            % get the data frame
            thisFrame = squeeze(perimeter.data(:,:,ii));
            
            % get the boundary points
            [Yc, Xc] = ind2sub(size(thisFrame),find(thisFrame));
            
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
                
                % loop over the parameters of the ellipse
                for jj=1:nEllipseParams
                    % Get the dataVector, restricted to the window range
                    dataVector=squeeze(loopVar_pLikelihoodMeanFromLastLoop(:,jj))';
                    dataVector=dataVector(rangeLowSignal:rangeHiSignal);
                    
                    % Build the precisionVector as the inverse of the measurement
                    % SD on each frame, scaled to range within the window from zero
                    % to unity. Thus, the noisiest measurement will not influence
                    % the prior.
                    precisionVector=squeeze(loopVar_pLikelihoodSDFromLastLoop(:,jj))';
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
                end % loop over ellipse params
                
                % The likelihood for this pass is the posterior from
                % the last loop
                pLikelihoodMeanTransparent = loopVar_pLikelihoodMeanFromLastLoop(ii,:);
                pLikelihoodSDTransparent = loopVar_pLikelihoodSDFromLastLoop(ii,:);
                
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
                
            end % not an empty frame
            
            % store results
            loopVar_pPosteriorMeanTransparent(ii,:) = pPosteriorMeanTransparent';
            loopVar_pPosteriorSDTransparent(ii,:) = pPosteriorSDTransparent';
            loopVar_finalFitError(ii) = fitError;
            
        end % parfor loop over frames
    end % loop over additional Bayes passes
end % check if we are doing additional Bayes passes


%% Clean up and save the fit results

% gather the loop vars into the ellipse structure
pupilData.pPriorMeanTransparent=loopVar_pPriorMeanTransparent;
pupilData.pPriorSDTransparent=loopVar_pPriorSDTransparent;
pupilData.pPosteriorMeanTransparent=loopVar_pPosteriorMeanTransparent;
pupilData.pPosteriorSDTransparent=loopVar_pPosteriorSDTransparent;
pupilData.fitError=loopVar_finalFitError';

% add a meta field with analysis details
pupilData.meta.smoothPupilParameters = p.Results;
pupilData.meta.smoothPupilParameters.coordinateSystem = 'intrinsicCoordinates(pixels)';

% save the pupilData
save(p.Results.pupilFileName,'pupilData')

% report completion of Bayesian analysis
if strcmp(p.Results.verbosity,'full')
    toc
    fprintf('\n');
end


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



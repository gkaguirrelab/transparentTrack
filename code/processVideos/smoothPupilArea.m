function [pupilData] = smoothPupilArea(perimeterFileName, pupilFileName, sceneGeometryFileName, varargin)
% smoothPupilArea(perimeterFileName, pupilFileName, sceneGeometryFileName)
%
% This routine implements a smoothing operation upon pupil area. The pupil
% area in the image plane is projected back to a model eye, using the
% parameters of scene geometry. This accounts for variation in the apparent
% area of the pupil at the image plane due to eye rotation. The
% reconstructed area of the pupil is then treated as the likelihood in an
% empirical Bayes approach that obtains the posterior estimate of the pupil
% area. A non-causal, exponentially weighted window of surrounding area
% values serves as a prior.
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
%   sceneGeometryFileName: full path to the .mat file that contains the
%    sceneGeometry structure for this acquisition.
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
%     3 : 2, this gives us an eccentricity UB threshold of ~0.75.
%  'exponentialTauParam' - The time constant (in video frames) of the
%     decaying exponential weighting function for pupil area.
%   'likelihoodErrorExponent' - The SD of the parameters estimated for each
%     frame are raised to this exponent, to either to weaken (>1) or
%     strengthen (<1) the influence of the current measure on the
%     posterior.
%   'whichLikelihoodSD' - The variance of the measured parameters for a
%     frame can be estimated using different methods. This setting controls
%     which of these is used to set the SD of the likelihood in the
%     calculation of the posterior. Valid values:
%       'pInitialFitHessianSD'
%       'pInitialFitSplitsSD'
%   'areaIdx' - The index of the ellipse parameters that holds pupil area.
%
% OUTPUTS:
%   pupilData: A structure with multiple fields corresponding to the
%     parameters, SDs, and errors of the initial and final ellipse fits.


%% Parse vargin for options passed here
p = inputParser; p.KeepUnmatched = true;

% Required
p.addRequired('perimeterFileName',@ischar);
p.addRequired('pupilFileName',@ischar);
p.addRequired('sceneGeometryFileName',@ischar);

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
p.addParameter('ellipseTransparentLB',[0, 0, 800, 0, -0.5*pi],@isnumeric);
p.addParameter('ellipseTransparentUB',[640,480,20000,0.75, 0.5*pi],@isnumeric);
p.addParameter('exponentialTauParam',5,@isnumeric);
p.addParameter('likelihoodErrorExponent',1.0,@isnumeric);
p.addParameter('whichLikelihoodSD','pInitialFitSplitsSD',@ischar);
p.addParameter('areaIdx',3,@isnumeric);


%% Parse and check the parameters
p.parse(perimeterFileName, pupilFileName, sceneGeometryFileName, varargin{:});

nEllipseParams=5; % 5 params in the transparent ellipse form

if length(p.Results.ellipseTransparentLB)~=nEllipseParams
    error('Wrong number of elements in ellipseTransparentLB');
end
if length(p.Results.ellipseTransparentUB)~=nEllipseParams
    error('Wrong number of elements in ellipseTransparentUB');
end
if sum(p.Results.ellipseTransparentUB>=p.Results.ellipseTransparentLB)~=nEllipseParams
    error('Lower bounds must be equal to or less than upper bounds');
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

% load the sceneGeometry structure
dataLoad=load(p.Results.sceneGeometryFileName);
sceneGeometry=dataLoad.sceneGeometry;
clear dataLoad

% Assemble a vector that contains the eye center of rotation and create a
% variable to hold the projectionModel string
eyeCenterOfRotation = [sceneGeometry.eyeCenter.X sceneGeometry.eyeCenter.Y sceneGeometry.eyeCenter.Z];
projectionModel = sceneGeometry.meta.projectionModel;

% Create a non-linear constraint for the ellipse fit
nonlinconst = @(transparentEllipseParams) constrainEllipseBySceneGeometry(...
    transparentEllipseParams, ...
    sceneGeometry);

% determine how many frames we will process
if p.Results.nFrames == Inf
    nFrames=size(perimeter.data,3);
else
    nFrames = p.Results.nFrames;
end

% convert the ellipse parameters in pupil data to eye azimuth, elevation,
% an pupil area
for ii = 1:nFrames
    [pupilAzi(ii), pupilEle(ii), pupilArea(ii)] = ...
        pupilProjection_inv( ...
        pupilData.pInitialFitTransparent(ii,:), ...
        eyeCenterOfRotation, ...
        projectionModel);
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

% Recast perimeter.data into a sliced cell array to reduce parfor
% broadcast overhead
frameCellArray = arrayfun(@(ii) {squeeze(perimeter.data(:,:,ii))},1:1:nFrames);

% Set-up other variables to be non-broadcast
verbosity = p.Results.verbosity;
areaIdx = p.Results.areaIdx;
whichLikelihoodSD = p.Results.whichLikelihoodSD;
likelihoodErrorExponent = p.Results.likelihoodErrorExponent;


%% Conduct empirical Bayes smoothing

% There are different measures available for the SD of the
% parameters of the initial fit. The parameter 'whichLikelihoodSD'
% controls which one of these is used for the likelihood. check that
% the requested SD measure is available
if ~isfield(pupilData,whichLikelihoodSD)
    error('The requested estimate of fit SD is not available in pupilData');
end

% Set up the decaying exponential weighting function. The relatively large
% window (8 times the time constant) is used to handle the case in
% which there is a stretch of missing data, in which case the long tails of
% the exponential can provide the prior.
window=ceil(max([p.Results.exponentialTauParam*8,8]));
windowSupport=1:1:window;
baseExpFunc=exp(-1/p.Results.exponentialTauParam*windowSupport);

% The weighting function is symmetric about the current time point. The
% current time point is excluded (set to nan)
exponentialWeights=[fliplr(baseExpFunc) NaN baseExpFunc];

% Alert the user
if strcmp(p.Results.verbosity,'full')
    tic
    fprintf(['Bayesian smoothing. Started ' char(datetime('now')) '\n']);
    fprintf('| 0                      50                   100%% |\n');
    fprintf('.\n');
end

% Loop through the frames
parfor (ii = 1:nFrames, nWorkers)
    
    % update progress
    if strcmp(verbosity,'full')
        if mod(ii,round(nFrames/50))==0
            fprintf('\b.\n');
        end
    end
    
    % initialize some variables
    pPosteriorMeanTransparent=NaN(1,nEllipseParams);
    pPosteriorFitError=NaN;
    
    % get the data frame
    thisFrame = frameCellArray{ii};
    
    % get the boundary points
    [Yc, Xc] = ind2sub(size(thisFrame),find(thisFrame));
    
    % if this frame has data, and the initial ellipse fit is not nan,
    % then proceed to calculate the posterior
    if ~isempty(Xc) &&  ~isempty(Yc) && sum(isnan(pupilData.pInitialFitTransparent(ii,:)))==0
        % Calculate the pupil area prior. The prior mean is given by the
        % surrounding area values, weighted by a decaying exponential in
        % time and the inverse of the standard deviation of each measure.
        % The prior standard deviation is weighted only by time.
        
        % A bit of fussing with the range here to handle the start and the
        % end of the data vector
        rangeLowSignal=max([ii-window,1]);
        rangeHiSignal=min([ii+window,nFrames]);
        restrictLowWindow= max([(ii-window-1)*-1,0]);
        restrictHiWindow = max([(nFrames-ii-window)*-1,0]);
        
        % Get the dataVector, restricted to the window range
        dataVector=pupilArea(rangeLowSignal:rangeHiSignal);
        
        % Build the precisionVector as the inverse of the measurement
        % SD on each frame, scaled to range within the window from zero
        % to unity. Thus, the noisiest measurement will not influence
        % the prior.
        precisionVector=squeeze(pupilData.(whichLikelihoodSD)(:,areaIdx))';
        precisionVector=precisionVector.^(-1);
        precisionVector=precisionVector(rangeLowSignal:rangeHiSignal);
        precisionVector=precisionVector-nanmin(precisionVector);
        precisionVector=precisionVector/nanmax(precisionVector);
        
        % The temporal weight vector is simply the exponential weights,
        % restricted to the available data widow
        temporalWeightVector = ...
            exponentialWeights(1+restrictLowWindow:end-restrictHiWindow);
        
        % Combine the precision and time weights, and calculate the
        % prior mean
        combinedWeightVector=precisionVector.*temporalWeightVector;
        priorPupilAreaMean = nansum(dataVector.*combinedWeightVector,2)./ ...
            nansum(combinedWeightVector(~isnan(dataVector)),2);
        
        % Obtain the standard deviation of the prior
        priorPupilAreaSD = nanstd(dataVector,temporalWeightVector);
        
        % Retrieve the initialFit for this frame
        likelihoodPupilAreaMean = pupilArea(ii);
        likelihoodPupilAreaSD = pupilData.(whichLikelihoodSD)(ii,areaIdx);
        
        % Raise the estimate of the SD from the initial fit to an
        % exponent. This is used to adjust the relative weighting of
        % the current frame realtive to the prior
        likelihoodPupilAreaSD = likelihoodPupilAreaSD .^ likelihoodErrorExponent;
        
        % Calculate the posterior values for the pupil fits, given the
        % likelihood and the prior
        posteriorPupilAreaMean = priorPupilAreaSD.^2.*likelihoodPupilAreaMean./(priorPupilAreaSD.^2+likelihoodPupilAreaSD.^2) + ...
            likelihoodPupilAreaSD.^2.*priorPupilAreaMean./(priorPupilAreaSD.^2+likelihoodPupilAreaSD.^2);
        
        % Calculate the SD of the posterior for kicks, but we don't use
        % this for anything.
        posteriorPupilAreaSD = sqrt((priorPupilAreaSD.^2.*likelihoodPupilAreaSD.^2) ./ ...
            (priorPupilAreaSD.^2+likelihoodPupilAreaSD.^2));
        
        % Convert the posterior pupil area in the eye back to area in the
        % image plane
        reconstructedTransparentEllipse = ...
            pupilProjection_fwd(pupilAzi(ii), pupilEle(ii), posteriorPupilAreaMean, eyeCenterOfRotation, projectionModel);
        
        % Pin the parameters are re-fit the ellipse to obtain the error
        % term
        lb_pin = pupilData.pInitialFitTransparent(ii,:);
        ub_pin = pupilData.pInitialFitTransparent(ii,:);
        lb_pin(areaIdx)=reconstructedTransparentEllipse(areaIdx);
        ub_pin(areaIdx)=reconstructedTransparentEllipse(areaIdx);
        [pPosteriorMeanTransparent, ~, pPosteriorFitError] = constrainedEllipseFit(Xc,Yc, lb_pin, ub_pin, nonlinconst);
        
    end % check if there are any perimeter points to fit
    
    % store results
    loopVar_pPosteriorMeanTransparent(ii,:) = pPosteriorMeanTransparent';
    loopVar_pPosteriorFitError(ii) = pPosteriorFitError;
    
end % loop over frames to calculate the posterior


%% Clean up and save the fit results

% gather the loop vars into the ellipse structure
pupilData.pPosteriorMeanTransparent=loopVar_pPosteriorMeanTransparent;
pupilData.pPosteriorFitError=loopVar_pPosteriorFitError';

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

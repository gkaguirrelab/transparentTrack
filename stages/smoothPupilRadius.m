function [pupilData] = smoothPupilRadius(perimeterFileName, pupilFileName, sceneGeometryFileName, varargin)
% Empirical Bayes smoothing of pupil radius in the scene
%
% Syntax:
%  [pupilData] = smoothPupilRadius(perimeterFileName, pupilFileName, sceneGeometryFileName)
%
% Description:
%   This routine implements a temporal smoothing operation upon the time
%   series of pupil radius values using an empirical Bayes approach. A
%   non-causal, exponentially weighted window of surrounding radius values
%   serves as an empirical prior. The posterior value of the radius is then
%   used as a constraint and the perimeter of the pupil in the image plane
%   is re-fit with the model eye.
%
% Notes:
%   Parallel pool - Controlled by the key/value pair 'useParallel'. The
%   routine should gracefully fall-back on serial processing if the
%   parallel pool is unavailable. Each worker requires ~8 GB of memory to
%   operate. It is important to keep total RAM usage below the physical
%   memory limit to prevent swapping and a dramatic slow down in
%   processing. To use the parallel pool with TbTb, provide the identity of
%   the repo name in the 'tbtbRepoName', which is then used to configure
%   the workers.
%
% Inputs:
%   perimeterFileName     - Full path to a .mat file that contains the
%                           perimeter data.
%   pupilFileName         - Full path to the .mat file that contains the
%                           pupil data to be smoothed. This file will be
%                           over-written by the output.
%   sceneGeometryFileName - Full path to the .mat file that contains the
%                           sceneGeometry to be used.
%
% Optional key/value pairs (display and I/O):
%  'verbose'              - Logical. Default false.
%
% Optional key/value pairs (flow control)
%  'nFrames'              - Analyze fewer than the total number of frames.
%  'useParallel'          - If set to true, use the Matlab parallel pool
%  'nWorkers'             - Specify the number of workers in the parallel
%                           pool. If undefined the default number will be
%                           used.
%
% Optional key/value pairs (environment)
%  'tbSnapshot'           - This should contain the output of the
%                           tbDeploymentSnapshot performed upon the result
%                           of the tbUse command. This documents the state
%                           of the system at the time of analysis.
%  'timestamp'            - AUTOMATIC; The current time and date
%  'username'             - AUTOMATIC; The user
%  'hostname'             - AUTOMATIC; The host
%
% Optional key/value pairs (fitting)
%  'glintFileName'        - Full path to a glint file. When available, the
%                           glint is used to constrain the eyePose that is
%                           found to fit the pupil perimeter.
%  'eyePoseLB'            - Lower bound on the eyePose
%  'eyePoseUB'            - Upper bound on the eyePose
%  'exponentialTauParam'  - The time constant (in video frames) of the
%                           decaying exponential weighting function for
%                           pupil radius.
%  'likelihoodErrorMultiplier' - The SD of the parameters estimated for
%                           each frame are computed as the product of the
%                           RMSE of the ellipse fits to the pupil perimeter
%                           points, the non-linear non-uniformity of the
%                           distribution of the points in space, and by
%                           this value. Typically set to ~4 to result in an
%                           SD of 1 when the fit of the points is good.
%                           Make this value larger to increase the
%                           influence of the prior, and smaller to increase
%                           the influence of the measurement.
%  'fitLabel'             - Identifies the field in pupilData that contains
%                           the ellipse fit params for which the search
%                           will be conducted.
%  'fixedPriorPupilRadius' - Scalar that provides the mean (in mm) of the
%                           expected radius of the pupil aperture during
%                           this acquisition. If set to empty (the default)
%                           the routine derives this values from the
%                           sceneConstrained fit results.
%  'relativeCameraPositionFileName' - Char. This is the full path to a
%                           relativeCameraPosition.mat file that provides
%                           the relative position of the camera at each
%                           video frame relative to the initial position of
%                           the camera.
%
% Outputs:
%   pupilData             - A structure with multiple fields corresponding
%                           to the parameters, SDs, and errors of the
%                           initial and final ellipse fits.
%

%% Parse vargin for options passed here
p = inputParser; p.KeepUnmatched = true; p.PartialMatching = false;

% Required
p.addRequired('perimeterFileName',@ischar);
p.addRequired('pupilFileName',@ischar);
p.addRequired('sceneGeometryFileName',@ischar);

% Optional display and I/O params
p.addParameter('verbose',false,@islogical);

% Optional flow control params
p.addParameter('nFrames',Inf,@isnumeric);
p.addParameter('useParallel',false,@islogical);
p.addParameter('nWorkers',[],@(x)(isempty(x) | isnumeric(x)));

% Optional environment parameters
p.addParameter('tbSnapshot',[],@(x)(isempty(x) | isstruct(x)));
p.addParameter('timestamp',char(datetime('now')),@ischar);
p.addParameter('hostname',char(java.lang.System.getProperty('user.name')),@ischar);
p.addParameter('username',char(java.net.InetAddress.getLocalHost.getHostName),@ischar);

% Optional fitting params
p.addParameter('glintFileName',[],@(x)(isempty(x) || ischar(x)));
p.addParameter('eyePoseLB',[-89,-89,0,0.1],@isnumeric);
p.addParameter('eyePoseUB',[89,89,0,5],@isnumeric);
p.addParameter('cameraTransBounds',[5; 5; 0],@isnumeric);
p.addParameter('exponentialTauParam',5,@isnumeric);
p.addParameter('likelihoodErrorMultiplier',1.0,@isnumeric);
p.addParameter('fixedPriorPupilRadius',[],@isnumeric);
p.addParameter('relativeCameraPositionFileName',[],@ischar);


%% Parse and check the parameters
p.parse(perimeterFileName, pupilFileName, sceneGeometryFileName, varargin{:});

nEllipseParams=5; % 5 params in the transparent ellipse form
nEyePoseParams=4; % 4 eyePose values (azimuth, elevation, torsion, radius)
radiusIdx = 4; % The 4th eyePose entry holds the radius value
nHeadTransParams=3; % [horizontal; vertical; depth]


%% Load and check data
% Load the pupil perimeter data
load(perimeterFileName,'perimeter');

% Load the pupil data
load(pupilFileName,'pupilData');

% Load the sceneGeometry file
% Need to do the inelegant loading approach to keep parloop happy
dataLoad = load(sceneGeometryFileName);
sceneGeometry = dataLoad.sceneGeometry;
clear dataLoad

% Load the glint file if passed
if ~isempty(p.Results.glintFileName)
    load(p.Results.glintFileName,'glintData');
    % Sometimes the glint file has fewer frames than the perimeter data. If
    % so, make sure that we discard any perimeter frames that do not have a
    % corresponding glint
    nGlintFrames = size(glintData.X,1);
    if nGlintFrames < size(perimeter.data,1)
        perimeter.data = perimeter.data(1:nGlintFrames);
    end
else
    glintData = [];
end


% Load the relativeCameraPosition file if passed and it exists
if ~isempty(p.Results.relativeCameraPositionFileName)
    if exist(p.Results.relativeCameraPositionFileName,'file')==2
        load(p.Results.relativeCameraPositionFileName,'relativeCameraPosition');
    else
        relativeCameraPosition=[];
    end
else
    relativeCameraPosition=[];
end

% Select the appropriate field of relativeCameraTransition, or synthesize
% one if not available
if ~isempty(relativeCameraPosition)
    cameraTransVec = ...
        relativeCameraPosition.sceneConstrained.values;
else
    cameraTransVec = zeros(nHeadTransParams,size(perimeter.data,1));
end

% determine how many frames we will process
if p.Results.nFrames == Inf
    nFrames=size(perimeter.data,1);
else
    nFrames = p.Results.nFrames;
end


%% Derive a fixed prior from a previous analysis of the data
% If not explicitly set, derive a prior across the entire acquisition for
% the mean and SD of the pupil size fom the sceneConstained results.
% Instead of the mean, we obtain the weighted median to avoid the influence
% of outlier values. The SD is set to something large (4 mm) so that the
% effect of this prior only appears in the near absence of other measures.
if isempty(p.Results.fixedPriorPupilRadius)
    fixedPriorPupilRadiusMean = medianw( ...
        pupilData.sceneConstrained.eyePoses.values(:,4), ...
        pupilData.sceneConstrained.ellipses.RMSE, 1 );
    fixedPriorPupilRadiusSD = 4;
else
    fixedPriorPupilRadiusMean = p.Results.fixedPriorPupilRadius;
    fixedPriorPupilRadiusSD = 4;
end


%% Set up the decaying exponential weighting function
% The relatively large window (10 times the time constant) is used to
% handle the case in which there is a stretch of missing data, in which
% case the long tails of the exponential can provide the prior.
window=ceil(max([p.Results.exponentialTauParam*10,10]));
windowSupport=1:1:window;
baseExpFunc=exp(-1/p.Results.exponentialTauParam*windowSupport);

% The weighting function is symmetric about the current time point. The
% current time point is excluded (set to nan)
exponentialWeights=[fliplr(baseExpFunc) NaN baseExpFunc];


%% Calculate likelhood SD across frames
% Obtain a measure for each frame of how completely the perimeter points
% define a full, 360 degrees around the pupil. This index of coverage of
% the pupil perimeter is distVals. If there is a perfectly uniform angular
% distribution of points in space around the pupil perimeter, then the
% distVals value will be zero. If there are perimeter points only at a
% single angular location around the pupil cirle, then the distVal will 1.

% The likelihood SD is based upon the RMSE of the fit of the elipse to the
% perimeter points for each frame
RMSE = pupilData.sceneConstrained.ellipses.RMSE';

% Define the bins over which the distribution of perimeter angles will be
% evaluated. 20 bins works pretty well.
nDivisions = 20;
histBins = linspace(-pi,pi,nDivisions);

% Anonymous function returns the linear non-uniformity of a set of values,
% ranging from 0 when perfectly uniform to 1 when completely non-uniform.
nonUniformity = @(x) (sum(abs(x/sum(x)-mean(x/sum(x))))/2)/(1-1/length(x));

% Loop over frames. Frames which have no perimeter points will be given a
% distVal of NaN.
for ii = 1:nFrames
    
    % Obtain the center of this fitted ellipse
    centerX = pupilData.sceneConstrained.ellipses.values(ii,1);
    centerY = pupilData.sceneConstrained.ellipses.values(ii,2);
    
    % Obtain the set of perimeter points
    Xp = perimeter.data{ii}.Xp;
    Yp = perimeter.data{ii}.Yp;
    
    % Calculate the deviation of the distribution of points from uniform
    linearNonUniformity(ii) = nonUniformity(histcounts(atan2(Yp-centerY,Xp-centerX),histBins));
end

% Subject the linearNonUniformity vector to a non-linear transformation.
% This has the effect of changing 0 --> 0.1, 0.8 --> 1, and values > 0.8
% --> infinity. This causes pupil perimeters with support at a single
% location (as opposed to fully around the perimeter) to have a markedly
% increased likelihood SD. Also, set InF values to something arbitrarily
% large.
distVals = (1./(1-sqrt(linearNonUniformity)))./10;
distVals(isinf(distVals)) = 1e20;

% The likelihood SD for each frame is the RMSE multiplied by the distVal
likelihoodPupilRadiusSDVector = distVals.*RMSE;

% Apply a multiplier that is used to adjust the relative weighting of the
% likelihood SD.
likelihoodPupilRadiusSDVector = likelihoodPupilRadiusSDVector .* p.Results.likelihoodErrorMultiplier;


%% Set up the parallel pool
if p.Results.useParallel
    nWorkers = startParpool( p.Results.nWorkers, p.Results.verbose );
else
    nWorkers=0;
end

% Recast perimeter.data into a sliced cell array to reduce parfor
% broadcast overhead
frameCellArray = perimeter.data(1:nFrames);
clear perimeter

% Set-up other variables to be non-broadcast
verbose = p.Results.verbose;
eyePoseLB = p.Results.eyePoseLB;
eyePoseUB = p.Results.eyePoseUB;
cameraTransBounds = p.Results.cameraTransBounds;


%% Perform the calculation across frames

% Alert the user
if p.Results.verbose
    tic
    fprintf(['Bayesian smoothing. Started ' char(datetime('now')) '\n']);
    fprintf('| 0                      50                   100%% |\n');
    fprintf('.\n');
end

% Store the warning state
warnState = warning();

% Loop through the frames
parfor (ii = 1:nFrames, nWorkers)
%for ii = 1:nFrames
    
    % update progress
    if verbose
        if mod(ii,round(nFrames/50))==0
            fprintf('\b.\n');
        end
    end
    
    % initialize some variables so that their use is transparent to the
    % parfor loop
    empiricalPriorPupilRadiusMean = NaN;
    empiricalPriorPupilRadiusSD = NaN;
    posteriorEllipseParams = NaN(1,nEllipseParams);
    posteriorEyePoseObjectiveError = NaN;
    posteriorEyePose = NaN(1,nEyePoseParams);
    posteriorCameraTrans = NaN(nHeadTransParams,1);
    posteriorPupilRadiusSD = NaN;
    uniformity = NaN;
    fitAtBound = false;
    
    % get the boundary points
    Xp = frameCellArray{ii}.Xp;
    Yp = frameCellArray{ii}.Yp;
    
    % Get the camera translation for this frame
    cameraTrans = cameraTransVec(:,ii);
    
    % if this frame has data, and eyePose radius is not nan, then proceed
    % to calculate the posterior
    if ~isempty(Xp) &&  ~isempty(Yp) && ~isnan(pupilData.sceneConstrained.eyePoses.values(ii,radiusIdx))
        % Calculate the pupil radius prior. The prior mean is given by the
        % surrounding radius values, weighted by a decaying exponential in
        % time and the inverse of the standard deviation of each measure.
        % The prior standard deviation is weighted only by time.
        
        % A bit of fussing with the range here to handle the start and the
        % end of the data vector
        rangeLowSignal=max([ii-window,1]);
        rangeHiSignal=min([ii+window,nFrames]);
        restrictLowWindow= max([(ii-window-1)*-1,0]);
        restrictHiWindow = max([(nFrames-ii-window)*-1,0]);
        
        % Get the dataVector, restricted to the window range
        dataVector=squeeze(pupilData.sceneConstrained.eyePoses.values(rangeLowSignal:rangeHiSignal,radiusIdx))';
        
        % The precisionVector is the inverse of the likelihood SD vector
        precisionVector = likelihoodPupilRadiusSDVector(rangeLowSignal:rangeHiSignal).^(-1);
        
        % Scale the precision vector within the window from zero to unity.
        % Thus, the noisiest measurement will not influence the prior.
        precisionVector=precisionVector-nanmin(precisionVector);
        precisionVector=precisionVector/nanmax(precisionVector);
        
        % The temporal weight vector is simply the exponential weights,
        % restricted to the available data widow
        temporalWeightVector = ...
            exponentialWeights(1+restrictLowWindow:end-restrictHiWindow);
        
        % Combine the precision and time weights, and calculate the
        % empirical prior pupil radius mean
        combinedWeightVector=precisionVector.*temporalWeightVector;
        empiricalPriorPupilRadiusMean = nansum(dataVector.*combinedWeightVector,2)./ ...
            nansum(combinedWeightVector(~isnan(dataVector)),2);
        
        % Obtain the standard deviation of the empirical prior, weighted
        % over time
        empiricalPriorPupilRadiusSD = nanstd(dataVector,temporalWeightVector);
        
        % Obtain the combined prior, which is the posterior of the fixed
        % and empirical priors
        combinedPriorPupilRadiusMean = fixedPriorPupilRadiusSD.^2.*empiricalPriorPupilRadiusMean./(fixedPriorPupilRadiusSD.^2+empiricalPriorPupilRadiusSD.^2) + ...
            empiricalPriorPupilRadiusSD.^2.*fixedPriorPupilRadiusMean./(fixedPriorPupilRadiusSD.^2+empiricalPriorPupilRadiusSD.^2);
        
        combinedPriorPupilRadiusSD = sqrt((fixedPriorPupilRadiusSD.^2.*empiricalPriorPupilRadiusSD.^2) ./ ...
            (fixedPriorPupilRadiusSD.^2+empiricalPriorPupilRadiusSD.^2));
        
        % Retrieve the initialFit for this frame
        likelihoodPupilRadiusMean = pupilData.sceneConstrained.eyePoses.values(ii,radiusIdx);
        likelihoodPupilRadiusSD = likelihoodPupilRadiusSDVector(ii);
        
        % Check if the likelihoodPupilRadiusSD is nan, in which case set it
        % to an arbitrarily large number so that the prior dictates the
        % posterior
        if isnan(likelihoodPupilRadiusSD)
            likelihoodPupilRadiusSD = 1e20;
        end
        
        % Calculate the posterior values for the pupil fits, given the
        % likelihood and the prior
        posteriorPupilRadius = combinedPriorPupilRadiusSD.^2.*likelihoodPupilRadiusMean./(combinedPriorPupilRadiusSD.^2+likelihoodPupilRadiusSD.^2) + ...
            likelihoodPupilRadiusSD.^2.*combinedPriorPupilRadiusMean./(combinedPriorPupilRadiusSD.^2+likelihoodPupilRadiusSD.^2);
        
        % Calculate the SD of the posterior of the pupil radius
        posteriorPupilRadiusSD = sqrt((combinedPriorPupilRadiusSD.^2.*likelihoodPupilRadiusSD.^2) ./ ...
            (combinedPriorPupilRadiusSD.^2+likelihoodPupilRadiusSD.^2));
        
        % It can be the case that the prior mean is nan, due to this frame
        % having a measurement, but all surrounding frames being bad.
        % Detect this case, and set the posterior to the likelihood.
        if isnan(empiricalPriorPupilRadiusMean)
            posteriorPupilRadius = likelihoodPupilRadiusMean;
            posteriorPupilRadiusSD = inf;
        end
        
        % Re-fit the ellipse with the radius constrained to the posterior
        % value. Pass the prior azimuth and elevation as x0.
        lb_pin = eyePoseLB;
        ub_pin = eyePoseUB;
        lb_pin(radiusIdx)=posteriorPupilRadius;
        ub_pin(radiusIdx)=posteriorPupilRadius;
        x0 = pupilData.sceneConstrained.eyePoses.values(ii,:);
        x0(radiusIdx)=posteriorPupilRadius;
        
        % If we have glintData, extract the glintCoord, and allow
        % non-zero bounds on the cameraTrans search. If no glint data,
        % then lock the cameraTransBounds to zero.
        if ~isempty(glintData)
            glintCoord = [glintData.X(ii,:), glintData.Y(ii,:)];
        else
            glintCoord = [];
        end
        
        if isempty(glintCoord) || any(isnan(glintCoord))
            thisFrameCameraTransBounds = [0; 0; 0];
        else
            thisFrameCameraTransBounds = cameraTransBounds;
        end
        
        % Turn off warnings that can arise when fitting bad frames
        warning('off','projectModelEye:rayTracingError');
        warning('off','projectModelEye:ellipseFitFailed');
        warning('off','gkaModelEye:pupilEllipseFit');
        
        % Perform the fit
        [posteriorEyePose, posteriorCameraTrans, posteriorEyePoseObjectiveError, posteriorEllipseParams, fitAtBound] = ...
            eyePoseEllipseFit(Xp, Yp, glintCoord, sceneGeometry, ...
            'cameraTransX0',cameraTrans,...
            'cameraTransBounds',thisFrameCameraTransBounds,...
            'eyePoseLB', lb_pin, 'eyePoseUB', ub_pin, 'eyePoseX0', x0);
        
        % Calculate the uniformity of the distribution of perimeter points
        % around the center of the fitted ellipse
        uniformity = 1-nonUniformity(histcounts(atan2(Yp-posteriorEllipseParams(2),Xp-posteriorEllipseParams(1)),histBins));
        
        % Restore the warning state
        warning(warnState);
        
    end % check if there are any perimeter points to fit
    
    % store results
    loopVar_empiricalPriorPupilRadiusMean(ii) = empiricalPriorPupilRadiusMean;
    loopVar_empiricalPriorPupilRadiusSD(ii) = empiricalPriorPupilRadiusSD;
    loopVar_posteriorEllipseParams(ii,:) = posteriorEllipseParams';
    loopVar_posteriorUniformity(ii) = uniformity;
    loopVar_posterioreyePosesObjectiveError(ii) = posteriorEyePoseObjectiveError;
    loopVar_fitAtBound(ii) = fitAtBound;
    loopVar_posteriorCameraTrans(ii,:) = posteriorCameraTrans;
    loopVar_posteriorEyePoses(ii,:) = posteriorEyePose;
    loopVar_posteriorPupilRadiusSD(ii) = posteriorPupilRadiusSD;
    
end % loop over frames to calculate the posterior

% report completion of Bayesian analysis
if p.Results.verbose
    toc
    fprintf('\n');
end

%% Clean up and save the fit results

% Clear out any prior results in the radiusSmoothed field
pupilData.radiusSmoothed = [];

% gather the loop vars into the ellipses field
pupilData.radiusSmoothed.ellipses.values=loopVar_posteriorEllipseParams;
pupilData.radiusSmoothed.ellipses.RMSE=loopVar_posterioreyePosesObjectiveError';
pupilData.radiusSmoothed.ellipses.uniformity=loopVar_posteriorUniformity';
pupilData.radiusSmoothed.ellipses.meta.ellipseForm = 'transparent';
pupilData.radiusSmoothed.ellipses.meta.labels = {'x','y','area','eccentricity','theta'};
pupilData.radiusSmoothed.ellipses.meta.units = {'pixels','pixels','squared pixels','non-linear eccentricity','rads'};
pupilData.radiusSmoothed.ellipses.meta.coordinateSystem = 'intrinsic image';
pupilData.radiusSmoothed.ellipses.meta.linearNonUniformity = linearNonUniformity;

% gather the loop vars into the eyePoses field
pupilData.radiusSmoothed.eyePoses.values=loopVar_posteriorEyePoses;
pupilData.radiusSmoothed.eyePoses.radiusSD=loopVar_posteriorPupilRadiusSD';
pupilData.radiusSmoothed.eyePoses.fitAtBound = loopVar_fitAtBound';
pupilData.radiusSmoothed.eyePoses.meta.labels = {'azimuth','elevation','torsion','pupil radius'};
pupilData.radiusSmoothed.eyePoses.meta.units = {'deg','deg','deg','mm'};
pupilData.radiusSmoothed.eyePoses.meta.coordinateSystem = 'head fixed (extrinsic)';
pupilData.radiusSmoothed.eyePoses.meta.empiricalPriorPupilRadiusMean = loopVar_empiricalPriorPupilRadiusMean;
pupilData.radiusSmoothed.eyePoses.meta.empiricalPriorPupilRadiusSD = loopVar_empiricalPriorPupilRadiusSD;
pupilData.radiusSmoothed.eyePoses.meta.likelihoodPupilRadiusSD = likelihoodPupilRadiusSDVector;
pupilData.radiusSmoothed.eyePoses.meta.fixedPriorPupilRadiusMean = fixedPriorPupilRadiusMean;
pupilData.radiusSmoothed.eyePoses.meta.fixedPriorPupilRadiusSD = fixedPriorPupilRadiusSD;

% add a meta field with analysis details
pupilData.radiusSmoothed.meta = p.Results;


% Update the relativeCameraPosition
relativeCameraPosition.radiusSmoothed.values = loopVar_posteriorCameraTrans';
relativeCameraPosition.radiusSmoothed.meta = p.Results;
relativeCameraPosition.currentField = 'radiusSmoothed';

% Store the identity of the most recently produced field of data
pupilData.currentField = 'radiusSmoothed';

% save the pupilData
save(pupilFileName,'pupilData')

% save the relativeCameraPosition
save(p.Results.relativeCameraPositionFileName,'relativeCameraPosition')

end % function

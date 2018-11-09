function [pupilData] = smoothPupilRadius(perimeterFileName, pupilFileName, sceneGeometryFileName, varargin)
% Empirical Bayes smoothing of pupil radius in the scene
%
% Syntax:
%  [pupilData] = smoothPupilRadius(perimeterFileName, pupilFileName, sceneGeometryFileName)
%
% Description:
%   This routine implements a temporal smoothing operation upon pupil
%   radius using an empirical Bayes approach. A non-causal, exponentially
%   weighted window of surrounding radius values serves as an empirical
%   prior. The posterior value of the radius is then used as a constraint
%   and the ellipse in the image plane is re-fit.
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
%  'eyePoseLB'            - Lower bound on the eyePose
%  'eyePoseUB'            - Upper bound on the eyePose
%  'exponentialTauParam'  - The time constant (in video frames) of the
%                           decaying exponential weighting function for
%                           pupil radius.
%  'likelihoodErrorMultiplier' - The SD of the parameters estimated for 
%                           each frame are multiplied by this value, to
%                           either to weaken (>1) or strengthen (<1) the
%                           influence of the current measure on the
%                           posterior.
%  'badFrameErrorThreshold' - Frames with RMSE fitting error above this
%                           threshold have their posterior values
%                           determined entirely by the prior. Additionally,
%                           these frames do not contribue to the prior.
%  'fitLabel'             - Identifies the field in pupilData that contains
%                           the ellipse fit params for which the search
%                           will be conducted.
%  'fixedPriorPupilRadius' - A 2x1 vector that provides the mean and SD (in 
%                           mm) of the expected radius of the pupil
%                           aperture during this acquisition. The defauklt
%                           values correspond to the pupil radius in a
%                           young adult in complete darkness.
%
% Outputs:
%   pupilData             - A structure with multiple fields corresponding
%                           to the parameters, SDs, and errors of the
%                           initial and final ellipse fits.
%

%% Parse vargin for options passed here
p = inputParser; p.KeepUnmatched = true;

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
p.addParameter('eyePoseLB',[-89,-89,0,0.1],@isnumeric);
p.addParameter('eyePoseUB',[89,89,0,5],@isnumeric);
p.addParameter('exponentialTauParam',3,@isnumeric);
p.addParameter('likelihoodErrorMultiplier',1.0,@isnumeric);
p.addParameter('badFrameErrorThreshold',2,@isnumeric);
p.addParameter('fitLabel','sceneConstrained',@ischar);
p.addParameter('fixedPriorPupilRadius',[3.5,1.5],@isnumeric);
p.addParameter('relativeCameraPositionFileName',[],@ischar);


%% Parse and check the parameters
p.parse(perimeterFileName, pupilFileName, sceneGeometryFileName, varargin{:});

nEllipseParams=5; % 5 params in the transparent ellipse form
nEyePoseParams=4; % 4 eyePose values (azimuth, elevation, torsion, radius) 
radiusIdx = 4; % The 4th eyePose entry holds the radius value

% Load the pupil perimeter data. It will be a structure variable
% "perimeter", with the fields .data and .meta
dataLoad=load(perimeterFileName);
perimeter=dataLoad.perimeter;
clear dataLoad

% Load the pupil data. It will be a structure variable "pupilData"
dataLoad=load(pupilFileName);
pupilData=dataLoad.pupilData;
clear dataLoad

% Load the sceneGeometry file
dataLoad=load(sceneGeometryFileName);
sceneGeometry=dataLoad.sceneGeometry;
clear dataLoad

% An earlier version of the code defined a non-zero iris thickness. We
% force this to zero here to speed computation
sceneGeometry.eye.iris.thickness=0;

% Load the relativeCameraPosition file if passed and it exists
if ~isempty(p.Results.relativeCameraPositionFileName)
    if exist(p.Results.pupilFileName, 'file')==2
        dataLoad=load(p.Results.relativeCameraPositionFileName);
        relativeCameraPosition=dataLoad.relativeCameraPosition;
        clear dataLoad
    else
        relativeCameraPosition=[];
    end
else
    relativeCameraPosition=[];
end

% determine how many frames we will process
if p.Results.nFrames == Inf
    nFrames=size(perimeter.data,1);
else
    nFrames = p.Results.nFrames;
end

% Check that the needed fields in the pupilData structure are present
if ~isfield(pupilData,(p.Results.fitLabel))
    error('The requested fit field is not available in pupilData');
end
if ~isfield(pupilData.(p.Results.fitLabel).ellipses,'RMSE')
    error('This fit field does not have the required subfield: ellipse.RMSE');
end


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
likelihoodErrorMultiplier = p.Results.likelihoodErrorMultiplier;
eyePoseLB = p.Results.eyePoseLB;
eyePoseUB = p.Results.eyePoseUB;
badFrameErrorThreshold = p.Results.badFrameErrorThreshold;
fitLabel = p.Results.fitLabel;
fixedPriorPupilRadiusMean = p.Results.fixedPriorPupilRadius(1);
fixedPriorPupilRadiusSD = p.Results.fixedPriorPupilRadius(2);

% Set up the decaying exponential weighting function. The relatively large
% window (10 times the time constant) is used to handle the case in which
% there is a stretch of missing data, in which case the long tails of the
% exponential can provide the prior.
window=ceil(max([p.Results.exponentialTauParam*10,10]));
windowSupport=1:1:window;
baseExpFunc=exp(-1/p.Results.exponentialTauParam*windowSupport);

% The weighting function is symmetric about the current time point. The
% current time point is excluded (set to nan)
exponentialWeights=[fliplr(baseExpFunc) NaN baseExpFunc];

% Obtain the RMSE of the fit of the elipse to the perimeter points for each
% frame
RMSE = pupilData.(p.Results.fitLabel).ellipses.RMSE';

% An older version of the code set a value of 1e12 for frames where the
% fitting failed. Just make these nans here.
RMSE(RMSE==1e12)=nan;

% Obtain a measure for each frame of how completely the perimeter points
% define the full circle of the pupil. If there is a perfectly uniform
% distribution of points around the pupil perimeter, then the distVals
% value will be zero. If there are perimeter points only at a single
% location around the pupil cirle, then the distVal will be ~4.35.
rmse = @(x) (sqrt(mean((x-mean(x)).^2)))/mean(x);
nDivisions = 20;
for ii = 1:nFrames
    distVals(ii) = rmse(histcounts(atan2(frameCellArray{ii}.Yp-pupilData.(p.Results.fitLabel).ellipses.values(ii,2),frameCellArray{ii}.Xp-pupilData.(p.Results.fitLabel).ellipses.values(ii,1)),linspace(-pi,pi,nDivisions)));
end
% Frames which have no perimeter points will return a distVal of zero. We
% set these to be nan.
distVals(distVals==0)=nan;

% The likelihood SD for each frame is the RMSE multiplied by the distVal
likelihoodPupilRadiusSDVector = distVals.*RMSE;


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
%for ii = 300:nFrames

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
    posteriorPupilRadiusSD = NaN;
    
    % get the boundary points
    Xp = frameCellArray{ii}.Xp;
    Yp = frameCellArray{ii}.Yp;
    
    % if this frame has data, and eyePose radius is not nan, then proceed
    % to calculate the posterior
    if ~isempty(Xp) &&  ~isempty(Yp) && ~isnan(pupilData.(fitLabel).eyePoses.values(ii,radiusIdx))
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
        dataVector=squeeze(pupilData.(fitLabel).eyePoses.values(rangeLowSignal:rangeHiSignal,radiusIdx))';
        
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
        likelihoodPupilRadiusMean = pupilData.(fitLabel).eyePoses.values(ii,radiusIdx);
        likelihoodPupilRadiusSD = likelihoodPupilRadiusSDVector(ii);
        
        % Apply a multiplier that is used to adjust the relative weighting
        % of the current frame relative to the prior
        likelihoodPupilRadiusSD = likelihoodPupilRadiusSD .* likelihoodErrorMultiplier;
                
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
        x0 = pupilData.(fitLabel).eyePoses.values(ii,:);
        x0(radiusIdx)=posteriorPupilRadius;
        
        % If a relativeCameraPosition is defined, update the sceneGeometry
        adjustedSceneGeometry = sceneGeometry;
        if ~isempty(relativeCameraPosition)
            cameraPosition = sceneGeometry.cameraPosition.translation;
            cameraPosition = cameraPosition + relativeCameraPosition(:,ii);
            adjustedSceneGeometry.cameraPosition.translation = cameraPosition;
        end        
        
        % Turn off warnings that can arise when fitting bad frames
        warning('off','pupilProjection_fwd:rayTracingError');
        warning('off','pupilProjection_fwd:ellipseFitFailed');
        
        % Perform the fit
        [posteriorEyePose, posteriorEyePoseObjectiveError, posteriorEllipseParams] = ...
            eyePoseEllipseFit(Xp, Yp, adjustedSceneGeometry, 'eyePoseLB', lb_pin, 'eyePoseUB', ub_pin, 'x0', x0, 'repeatSearchThresh', badFrameErrorThreshold);
        
        % Restore the warning state
        warning(warnState);
        
    end % check if there are any perimeter points to fit
    
    % store results
    loopVar_empiricalPriorPupilRadiusMean(ii) = empiricalPriorPupilRadiusMean;
    loopVar_empiricalPriorPupilRadiusSD(ii) = empiricalPriorPupilRadiusSD;
    loopVar_posteriorEllipseParams(ii,:) = posteriorEllipseParams';
    loopVar_posterioreyePosesObjectiveError(ii) = posteriorEyePoseObjectiveError;
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
pupilData.radiusSmoothed.ellipses.meta.ellipseForm = 'transparent';
pupilData.radiusSmoothed.ellipses.meta.labels = {'x','y','area','eccentricity','theta'};
pupilData.radiusSmoothed.ellipses.meta.units = {'pixels','pixels','squared pixels','non-linear eccentricity','rads'};
pupilData.radiusSmoothed.ellipses.meta.coordinateSystem = 'intrinsic image';

% gather the loop vars into the eyePoses field
pupilData.radiusSmoothed.eyePoses.values=loopVar_posteriorEyePoses;
pupilData.radiusSmoothed.eyePoses.radiusSD=loopVar_posteriorPupilRadiusSD';
pupilData.radiusSmoothed.eyePoses.meta.labels = {'azimuth','elevation','torsion','pupil radius'};
pupilData.radiusSmoothed.eyePoses.meta.units = {'deg','deg','deg','mm'};
pupilData.radiusSmoothed.eyePoses.meta.coordinateSystem = 'head fixed (extrinsic)';
pupilData.radiusSmoothed.eyePoses.meta.empiricalPriorPupilRadiusMean = loopVar_empiricalPriorPupilRadiusMean;
pupilData.radiusSmoothed.eyePoses.meta.empiricalPriorPupilRadiusSD = loopVar_empiricalPriorPupilRadiusSD;

% add a meta field with analysis details
pupilData.radiusSmoothed.meta = p.Results;

% save the pupilData
save(p.Results.pupilFileName,'pupilData')

end % function

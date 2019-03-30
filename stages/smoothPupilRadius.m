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
%  'fitLabel'             - Identifies the field in pupilData that contains
%                           the ellipse fit params for which the search
%                           will be conducted.
%  'fixedPriorPupilRadius' - A 2x1 vector that provides the mean and SD (in 
%                           mm) of the expected radius of the pupil
%                           aperture during this acquisition. The default
%                           values correspond to the pupil radius in a
%                           young adult in complete darkness.
%  'adjustedCameraPositionTranslation' - 3x1 vector that provides position
%                           of the camera relative to the origin of the
%                           world coordinate system (which is the anterior
%                           surface of the cornea in primary gaze). This
%                           value is used to update the sceneGeometry file
%                           to account for head movement that has taken
%                           place between the sceneGeometry acquisition and
%                           the acquisition undergoing analysis. This
%                           updated camera position should reflect the
%                           camera position at the start of the current
%                           acquisition.
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
p.addParameter('likelihoodErrorMultiplier',4.0,@isnumeric);
p.addParameter('fitLabel','sceneConstrained',@ischar);
p.addParameter('fixedPriorPupilRadius',[3.5,1.5],@isnumeric);
p.addParameter('adjustedCameraPositionTranslation',[],@isnumeric);
p.addParameter('relativeCameraPositionFileName',[],@ischar);


%% Parse and check the parameters
p.parse(perimeterFileName, pupilFileName, sceneGeometryFileName, varargin{:});

nEllipseParams=5; % 5 params in the transparent ellipse form
nEyePoseParams=4; % 4 eyePose values (azimuth, elevation, torsion, radius) 
radiusIdx = 4; % The 4th eyePose entry holds the radius value


%% Load and check data
% Load the pupil perimeter data
dataLoad=load(perimeterFileName);
perimeter=dataLoad.perimeter;
clear dataLoad

% Load the pupil data
dataLoad=load(pupilFileName);
pupilData=dataLoad.pupilData;
clear dataLoad

% Load the sceneGeometry file
dataLoad=load(sceneGeometryFileName);
sceneGeometry=dataLoad.sceneGeometry;
clear dataLoad

%% HACK
% An earlier version of the code defined a non-zero iris thickness. We
% force this to zero here to speed computation. This line should ultimately
% be removed when Geoff is done processing the session 1 Connectome data.
sceneGeometry.eye.iris.thickness=0;

% If an adjustedCameraPositionTranslation value has been passed, update
% this field of the sceneGeometry
if ~isempty(p.Results.adjustedCameraPositionTranslation)
    sceneGeometry.cameraPosition.translation = p.Results.adjustedCameraPositionTranslation;
end

% Load the relativeCameraPosition file if passed and it exists
if ~isempty(p.Results.relativeCameraPositionFileName)
    if exist(p.Results.relativeCameraPositionFileName, 'file')==2
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
RMSE = pupilData.(p.Results.fitLabel).ellipses.RMSE';

%% HACK
% An older version of the code set a value of 1e12 for frames where the
% fitting failed. Just make these nans here.
% Delete this line once Geoff has processed the session 1 TOME data.
RMSE(RMSE==1e12)=nan;

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
    centerX = pupilData.(p.Results.fitLabel).ellipses.values(ii,1);
    centerY = pupilData.(p.Results.fitLabel).ellipses.values(ii,2);

    % Obtain the set of perimeter points
    Xp = perimeter.data{ii}.Xp;
    Yp = perimeter.data{ii}.Xp;
    
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
fitLabel = p.Results.fitLabel;
fixedPriorPupilRadiusMean = p.Results.fixedPriorPupilRadius(1);
fixedPriorPupilRadiusSD = p.Results.fixedPriorPupilRadius(2);


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
    posteriorPupilRadiusSD = NaN;
    uniformity = NaN;
    fitAtBound = false;
    
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
            cameraPosition = cameraPosition - relativeCameraPosition.values(:,ii);
            adjustedSceneGeometry.cameraPosition.translation = cameraPosition;
        end
        
        % Turn off warnings that can arise when fitting bad frames
        warning('off','pupilProjection_fwd:rayTracingError');
        warning('off','pupilProjection_fwd:ellipseFitFailed');
        
        % Perform the fit
        [posteriorEyePose, posteriorEyePoseObjectiveError, posteriorEllipseParams, fitAtBound] = ...
            eyePoseEllipseFit(Xp, Yp, adjustedSceneGeometry, 'eyePoseLB', lb_pin, 'eyePoseUB', ub_pin, 'x0', x0);
        
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


% add a meta field with analysis details
pupilData.radiusSmoothed.meta = p.Results;

% save the pupilData
save(p.Results.pupilFileName,'pupilData')

end % function

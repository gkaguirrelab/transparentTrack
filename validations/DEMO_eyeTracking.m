% DEMO_eyeTracking
%
% Demonstrate the entire eyetracking analysis pipeline.
%
% A local sandbox folder named 'eyeTrackingDEMO' will be created on the
% desktop to replicate the dropbox environment of the real routine. Files
% will be downloaded from figshare and placed in the sandbox (about 7 GB).
%
% Make sure your machine is configured to work with ToolboxToolbox.
%
% Run-time on an average computer is about 5 minutes for 500 frames.
% Set nFrames to 'Inf' to process the entire video, which will take longer.
%
% Usage examples
% ==============
%
% DEMO_eyeTracking
%


%% hard coded parameters
nFrames = Inf; % number of frames to process (set to Inf to do all)
verbose = true; % Set to none to make the demo silent
TbTbToolboxName = 'transparentTrack';


%% set paths and make directories
% create test sandbox on desktop
sandboxDir = '~/Desktop/eyeTrackingDEMO';
if ~exist(sandboxDir,'dir')
    mkdir(sandboxDir)
end

% define path parameters
pathParams.dataSourceDirRoot = fullfile(sandboxDir,'TOME_data');
pathParams.dataOutputDirRoot = fullfile(sandboxDir,'TOME_processing');
pathParams.projectSubfolder = 'session2_spatialStimuli';
pathParams.eyeTrackingDir = 'EyeTracking';
pathParams.subjectID = 'TOME_3020';
pathParams.sessionDate = '050517';
pathParams.runName = 'GazeCal01';


%% TbTb configuration
% We will suppress the verbose output, but detect if there are deploy
% errors and if so stop execution
tbConfigResult=tbUse(TbTbToolboxName,'reset','full','verbose',false);
if sum(cellfun(@sum,extractfield(tbConfigResult, 'isOk')))~=length(tbConfigResult)
    error('There was a tb deploy error. Check the contents of tbConfigResult');
end
% We save a deployment snapshot. This variable is passed to the analysis
% pipeline and then saved with every output file, thereby documenting the
% system and software configuration at the time of execution.
tbSnapshot=tbDeploymentSnapshot(tbConfigResult,'verbose',false);
clear tbConfigResult


%% Prepare paths and directories
% define full paths for input and output
pathParams.dataSourceDirFull = fullfile(pathParams.dataSourceDirRoot, pathParams.projectSubfolder, ...
    pathParams.subjectID, pathParams.sessionDate, pathParams.eyeTrackingDir);
pathParams.dataOutputDirFull = fullfile(pathParams.dataOutputDirRoot, pathParams.projectSubfolder, ...
    pathParams.subjectID, pathParams.sessionDate, pathParams.eyeTrackingDir);

% Download the data if it is not already there
demoPackage = fullfile(sandboxDir,'eyeTrackingDEMO.zip');
if ~exist (demoPackage,'file')
    url = 'https://ndownloader.figshare.com/files/9355459?private_link=011191afe46841d2c2f5';
    system (['curl -L ' sprintf(url) ' > ' sprintf(demoPackage)])
    currentDir = pwd;
    cd (sandboxDir)
    unzip(demoPackage)
    cd (currentDir)
end


%% Prepare analysis parameters
% We can conceptually break the analysis into two components:
% - Initital processing to the stage of fitting an ellipse to the pupil
%   perimeter
% - Definition of the scene geometry and re-fitting of the pupil
%   perimeter using scene constraints and empirical Bayesian temporal 
%   smoothing.


%% Initial processing -- deinterlace to minimally constrained ellipse fitting
% These parameters are used in the function findPupilPerimeter; definition
% of the parameters may be found in the header comments for that routine.
% To select good parameters for your video, use the interactive routine 
% estimateSceneParamsGUI.m, which is found in the Utilities directory
pupilFrameMask = [64 109 75 183];
glintFrameMask = [157 148 173 192];
pupilRange = [34 51];
pupilCircleThresh = 0.0179;
pupilGammaCorrection = 0.75;

%% Run the analysis pipeline
% This routine will produce the initial ellipse fit and a fit video. It
% includes the stages:
%   deinterlaceVideo
%   findGlint
%   findPupilPerimeter
%   makeControlFile
%   applyControlFile
%   fitPupilPerimeter
%   makeFitVideo
%
runVideoPipeline( pathParams, ...
    'nFrames',nFrames,'verbose', verbose, 'tbSnapshot',tbSnapshot, 'useParallel',true, ...
    'pupilFrameMask', [64 109 75 183], 'glintFrameMask', [157 148 173 192], ...
    'pupilRange', [34 51], 'pupilCircleThresh', 0.0179, 'pupilGammaCorrection', 0.75, ...
    'overwriteControlFile', true, 'catchErrors', false,...
    'skipStageByNumber',[],'makeFitVideoByNumber',[6]);


%% Secondary processing -- definition and fitting with scene geometry
% Improved fitting can be obtained by defining properties of the eye, the
% camera, and their relationship. Parameters that define the camera and the
% eye described in greater detail in the routine createSceneGeometry which
% is found in the gkaModelEye (https://github.com/gkaguirrelab/gkaModelEye)
% code repository.

% Define camera parameters. These were obtained by an empirical measurement
% (camera resectioning) of the IR camera used to record the demo data. Use
% the matlab routine cameraCalibrator: 
%	https://www.mathworks.com/help/vision/ug/single-camera-calibrator-app.html
intrinsicCameraMatrix = [2627.0 0 338.1; 0 2628.1 246.2; 0 0 1];
radialDistortionVector = [-0.3517 3.5353];
sensorResolution = [640 480];
spectralDomain = 'nir';

% Define properties of the eye of the subject. The "maxIrisDiamPixels" is
% the largest horizontal diameter of the iris (colored portion) of the eye
% that is seen in the video. This valye is used to estimate the distance of
% the camera from the eye. The spherical ametropia is the refractive error
% correction (in diopters) of the eye of the subject. A negative value is
% the correction for a myopic (near-sighted) person. The axial length would
% be obtained from an ophthalmologic device such as the IOL Master, as
% would the measuredCornealCurvature. If any of these last three values are
% not available, leave these parameters undefined.
eyeLaterality = 'right';
maxIrisDiamPixels = 270;
sphericalAmetropia = -1.5;
axialLength = 25.35;
measuredCornealCurvature = [41.36,41.67,25];

% Estimate camera distance from iris diameter in pixels. Because biological
% variation in the size of the visible iris is known, we can use the
% observed maximum diameter of the iris in pixels to obtain a guess as to
% the distance of the eye from the camera.
sceneGeometry = createSceneGeometry(...
    'radialDistortionVector',radialDistortionVector, ...
    'intrinsicCameraMatrix',intrinsicCameraMatrix);
[cameraDepthMean, cameraDepthSD] = depthFromIrisDiameter( sceneGeometry, maxIrisDiamPixels );

% Assemble the scene parameter bounds. These are in the order of:
%   [torsion; x; y; z; eyeRotationScalarJoint; eyeRotationScalerDifferential]
% where torsion specifies the torsion of the camera with respect to the eye
% in degrees, [x y z] is the translation of the camera w.r.t. the eye in
% mm, and the eyeRotationScalar variables are multipliers that act upon the
% centers of rotation estimated for the eye.
% If the eye is markedly off-center in the image, then the translation
% bounds should be increased.
sceneParamsLB = [-5; -5; -5; cameraDepthMean-cameraDepthSD*2; 0.75; 0.9];
sceneParamsLBp = [-3; -2; -2; cameraDepthMean-cameraDepthSD*1; 0.85; 0.95];
sceneParamsUBp = [3; 2; 2; cameraDepthMean+cameraDepthSD*1; 1.15; 1.05];
sceneParamsUB = [5; 5; 5; cameraDepthMean+cameraDepthSD*2; 1.25; 1.1];

% The estimation of scene geometry is greatly aided by having the subject
% fixate targets at known visual angle positions. The routine
% estimateSceneParams acceps a list of frames of the video during which the
% pupil is well visualized and a list of target positions (X,Y) in units of
% degrees corresponding to each of those video frames.
fixationTargetArray = [ -7, 0, 7, 0, 0, 7, -7, 7, -7 ; 7, -7, 0, 0, 7, 7, 0, -7, -7];
ellipseArrayList = [ 732, 896, 1023, 1167, 1261, 1383, 1542, 1646, 1808 ];

% Run the video pipeline from the stage of estimation of scene geometry
% through to the end. It includes the stages:
%   estimateSceneParams
%   fitPupilPerimeter (now with scene geometry constraints)
%   smoothPupilRadius
%   makeFitVideo
%
runVideoPipeline( pathParams, ...
    'nFrames',nFrames,'verbose', verbose, 'tbSnapshot',tbSnapshot, 'useParallel',true, ...
    'intrinsicCameraMatrix',intrinsicCameraMatrix, ...
    'radialDistortionVector',radialDistortionVector, ...
    'sensorResolution',sensorResolution,...
    'spectralDomain',spectralDomain, ...
    'eyeLaterality',eyeLaterality,'sphericalAmetropia',sphericalAmetropia,...
    'axialLength',axialLength,'measuredCornealCurvature',measuredCornealCurvature,...
    'sceneParamsLB',sceneParamsLB,'sceneParamsUB',sceneParamsUB,...
    'sceneParamsLBp',sceneParamsLBp,'sceneParamsUBp',sceneParamsUBp,...
    'catchErrors', false,...
    'fixationTargetArray',fixationTargetArray,'ellipseArrayList', ellipseArrayList, ...
    'skipStageByNumber',[1:6]);




%% Plot some fits
pupilFileName = fullfile(pathParams.dataOutputDirFull,[pathParams.runName '_pupil.mat']);
dataLoad = load(pupilFileName);
pupilData = dataLoad.pupilData;
clear dataLoad

temporalSupport = 0:1/60.:(size(pupilData.sceneConstrained.ellipses.values,1)-1)/60; % seconds
temporalSupport = temporalSupport / 60; % minutes

% Points with a narrow posterior distribution of the fit to pupil radius
goodIdx = pupilData.radiusSmoothed.eyePoses.radiusSD < 0.01;

% Make a plot of pupil area, both on the image plane and on the eye
figure
subplot(2,1,1)
plot(temporalSupport,pupilData.initial.ellipses.values(:,3),'-k','LineWidth',2);
hold on
plot(temporalSupport,pupilData.sceneConstrained.ellipses.values(:,3),'-b');
plot(temporalSupport,pupilData.sceneConstrained.ellipses.values(:,3)-pupilData.sceneConstrained.ellipses.splitsSD(:,3),'-','Color',[0 0 0.7])
plot(temporalSupport,pupilData.sceneConstrained.ellipses.values(:,3)+pupilData.sceneConstrained.ellipses.splitsSD(:,3),'-','Color',[0 0 0.7])
plot(temporalSupport(goodIdx),pupilData.radiusSmoothed.ellipses.values(goodIdx,3),'.r','LineWidth',2)
xlim([0 max(temporalSupport)]);
xlabel('time [mins]');
ylabel('pupil area [pixels in plane]');
ylim([5000 25000]);
hold off

subplot(2,1,2)
plot(temporalSupport,pupilData.sceneConstrained.eyePoses.values(:,4),'-k','LineWidth',2);
hold on
plot(temporalSupport,pupilData.sceneConstrained.eyePoses.values(:,4),'-b');
plot(temporalSupport,pupilData.sceneConstrained.eyePoses.values(:,4)-pupilData.sceneConstrained.eyePoses.splitsSD(:,4),'-','Color',[0 0 0.7])
plot(temporalSupport,pupilData.sceneConstrained.eyePoses.values(:,4)+pupilData.sceneConstrained.eyePoses.splitsSD(:,4),'-','Color',[0 0 0.7])
plot(temporalSupport(goodIdx),pupilData.radiusSmoothed.eyePoses.values(goodIdx,4),'.r','LineWidth',2)
xlim([0 max(temporalSupport)]);
xlabel('time [mins]');
ylabel('pupil radius [mm on eye]');
ylim([0 4]);
hold off

% Make a plot of X and Y eye pupil position on the image plane
figure
subplot(2,1,1)
plot(temporalSupport,pupilData.initial.ellipses.values(:,1),'-k','LineWidth',2);
hold on
plot(temporalSupport,pupilData.sceneConstrained.ellipses.values(:,1),'-b');
plot(temporalSupport,pupilData.sceneConstrained.ellipses.values(:,1)-pupilData.sceneConstrained.ellipses.splitsSD(:,1),'-','Color',[0 0 0.7])
plot(temporalSupport,pupilData.sceneConstrained.ellipses.values(:,1)+pupilData.sceneConstrained.ellipses.splitsSD(:,1),'-','Color',[0 0 0.7])
plot(temporalSupport(goodIdx),pupilData.radiusSmoothed.ellipses.values(goodIdx,1),'.r','LineWidth',2)
xlim([0 max(temporalSupport)]);
xlabel('time [mins]');
ylabel('X position [pixels]');
ylim([0 500]);
hold off

subplot(2,1,2)
plot(temporalSupport,pupilData.initial.ellipses.values(:,2),'-k','LineWidth',2);
hold on
plot(temporalSupport,pupilData.sceneConstrained.ellipses.values(:,2),'-b');
plot(temporalSupport,pupilData.sceneConstrained.ellipses.values(:,2)-pupilData.sceneConstrained.ellipses.splitsSD(:,2),'-','Color',[0 0 0.7])
plot(temporalSupport,pupilData.sceneConstrained.ellipses.values(:,2)+pupilData.sceneConstrained.ellipses.splitsSD(:,2),'-','Color',[0 0 0.7])
plot(temporalSupport(goodIdx),pupilData.radiusSmoothed.ellipses.values(goodIdx,2),'.r','LineWidth',2)

xlim([0 max(temporalSupport)]);
xlabel('time [mins]');
ylabel('Y position [pixels]');
ylim([0 500]);
hold off
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

%% set paths and make directories
% create test sandbox on desktop
sandboxDir = '~/Desktop/eyeTrackingDEMO';
if ~exist(sandboxDir,'dir')
    mkdir(sandboxDir)
end

%% hard coded parameters
nFrames = Inf; % number of frames to process (set to Inf to do all)
verbosity = 'full'; % Set to none to make the demo silent
TbTbToolboxName = 'transparentTrack';

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

% Define intrinsic camera parameters
intrinsicCameraMatrix = [2627.0 0 338.1; 0 2628.1 246.2; 0 0 1];
radialDistortionVector = [-0.3517 3.5353];

% Define subject parameters
eyeLaterality = 'right';
axialLength = 25.35;
sphericalAmetropia = -1.5;

% Estimate camera distance from iris diameter in pixels
observedIrisDiamPixels = 267;
sceneGeometry = createSceneGeometry(...
    'radialDistortionVector',radialDistortionVector, ...
    'intrinsicCameraMatrix',intrinsicCameraMatrix, ...
    'constraintTolerance',0.02);
[cameraDepthMean, cameraDepthSD] = depthFromIrisDiameter( sceneGeometry, observedIrisDiamPixels );

% Assemble the scene parameter bounds
sceneParamsLB = [-5; -5; -5; cameraDepthMean-cameraDepthSD*2; 0.75; 0.9];
sceneParamsUB = [5; 5; 5; cameraDepthMean+cameraDepthSD*2; 1.25; 1.1];
sceneParamsLBp = [-3; -2; -2; cameraDepthMean-cameraDepthSD*1; 0.85; 0.95];
sceneParamsUBp = [3; 2; 2; cameraDepthMean+cameraDepthSD*1; 1.15; 1.05];

% Run the analysis pipeline
runVideoPipeline( pathParams, ...
    'nFrames',nFrames,'verbosity', verbosity, 'tbSnapshot',tbSnapshot, 'useParallel',true, ...
    'pupilRange', [40 200], 'pupilCircleThresh', 0.04, 'pupilGammaCorrection', 1.5, ...
    'intrinsicCameraMatrix',intrinsicCameraMatrix, ...
    'radialDistortionVector',radialDistortionVector, ...
    'eyeLaterality',eyeLaterality,'axialLength',axialLength,'sphericalAmetropia',sphericalAmetropia,...
    'sceneParamsLB',sceneParamsLB,'sceneParamsUB',sceneParamsUB,...
    'sceneParamsLBp',sceneParamsLBp,'sceneParamsUBp',sceneParamsUBp,...
    'overwriteControlFile', true, 'catchErrors', false,...
    'skipStageByNumber',[1:5],'makeFitVideoByNumber',[6 8]);


%% Plot some fits
pupilFileName = fullfile(pathParams.dataOutputDirFull,[pathParams.runName '_pupil.mat']);
dataLoad = load(pupilFileName);
pupilData = dataLoad.pupilData;
clear dataLoad

temporalSupport = 0:1/60.:(size(pupilData.sceneConstrained.ellipses.values,1)-1)/60; % seconds
temporalSupport = temporalSupport / 60; % minutes

% Make a plot of pupil area, both on the image plane and on the eye
figure
subplot(2,1,1)
plot(temporalSupport,pupilData.initial.ellipses.values(:,3),'-k','LineWidth',2);
hold on
plot(temporalSupport,pupilData.sceneConstrained.ellipses.values(:,3),'-b');
plot(temporalSupport,pupilData.sceneConstrained.ellipses.values(:,3)-pupilData.sceneConstrained.ellipses.splitsSD(:,3),'-','Color',[0 0 0.7])
plot(temporalSupport,pupilData.sceneConstrained.ellipses.values(:,3)+pupilData.sceneConstrained.ellipses.splitsSD(:,3),'-','Color',[0 0 0.7])
plot(temporalSupport,pupilData.radiusSmoothed.ellipses.values(:,3),'-r','LineWidth',2)
xlim([0 max(temporalSupport)]);
xlabel('time [mins]');
ylabel('pupil area [pixels in plane]');
hold off

subplot(2,1,2)
plot(temporalSupport,pupilData.sceneConstrained.eyePoses.values(:,4),'-k','LineWidth',2);
hold on
plot(temporalSupport,pupilData.sceneConstrained.eyePoses.values(:,4),'-b');
plot(temporalSupport,pupilData.sceneConstrained.eyePoses.values(:,4)-pupilData.sceneConstrained.eyePoses.splitsSD(:,4),'-','Color',[0 0 0.7])
plot(temporalSupport,pupilData.sceneConstrained.eyePoses.values(:,4)+pupilData.sceneConstrained.eyePoses.splitsSD(:,4),'-','Color',[0 0 0.7])
plot(temporalSupport,pupilData.radiusSmoothed.eyePoses.values(:,4),'-r','LineWidth',2)
xlim([0 max(temporalSupport)]);
xlabel('time [mins]');
ylabel('pupil radius [mm on eye]');
hold off

% Make a plot of X and Y eye pupil position on the image plane
figure
subplot(2,1,1)
plot(temporalSupport,pupilData.initial.ellipses.values(:,1),'-k','LineWidth',2);
hold on
plot(temporalSupport,pupilData.sceneConstrained.ellipses.values(:,1),'-b');
plot(temporalSupport,pupilData.sceneConstrained.ellipses.values(:,1)-pupilData.sceneConstrained.ellipses.splitsSD(:,1),'-','Color',[0 0 0.7])
plot(temporalSupport,pupilData.sceneConstrained.ellipses.values(:,1)+pupilData.sceneConstrained.ellipses.splitsSD(:,1),'-','Color',[0 0 0.7])
plot(temporalSupport,pupilData.radiusSmoothed.ellipses.values(:,1),'-r','LineWidth',2)
xlim([0 max(temporalSupport)]);
xlabel('time [mins]');
ylabel('X position [pixels]');
hold off

subplot(2,1,2)
plot(temporalSupport,pupilData.initial.ellipses.values(:,2),'-k','LineWidth',2);
hold on
plot(temporalSupport,pupilData.sceneConstrained.ellipses.values(:,2),'-b');
plot(temporalSupport,pupilData.sceneConstrained.ellipses.values(:,2)-pupilData.sceneConstrained.ellipses.splitsSD(:,2),'-','Color',[0 0 0.7])
plot(temporalSupport,pupilData.sceneConstrained.ellipses.values(:,2)+pupilData.sceneConstrained.ellipses.splitsSD(:,2),'-','Color',[0 0 0.7])
plot(temporalSupport,pupilData.radiusSmoothed.ellipses.values(:,2),'-r','LineWidth',2)

xlim([0 max(temporalSupport)]);
xlabel('time [mins]');
ylabel('Y position [pixels]');
hold off
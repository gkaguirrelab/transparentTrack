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


%% Perform the entire analysis with one call
runVideoPipeline( pathParams, ...
    'nFrames',nFrames,'verbosity', verbosity, 'tbSnapshot',tbSnapshot, 'useParallel',true, ...
    'pupilRange', [40 200], 'pupilCircleThresh', 0.04, 'pupilGammaCorrection', 1.5, ...
    'overwriteControlFile', true, 'catchErrors', false);


%% Plot some fits
pupilFileName = fullfile(pathParams.dataOutputDirFull,[pathParams.runName '_pupil.mat']);
dataLoad = load(pupilFileName);
pupilData = dataLoad.pupilData;
clear dataLoad

temporalSupport = 0:1/60.:(size(pupilData.ellipseParamsSceneConstrained_mean,1)-1)/60; % seconds
temporalSupport = temporalSupport / 60; % minutes

% Make a plot of pupil area, both on the image plane and on the eye
figure
subplot(2,1,1)
plot(temporalSupport,pupilData.ellipseParamsUnconstrained_mean(:,3),'-k','LineWidth',2);
hold on
plot(temporalSupport,pupilData.ellipseParamsSceneConstrained_mean(:,3),'-b');
plot(temporalSupport,pupilData.ellipseParamsSceneConstrained_mean(:,3)-pupilData.ellipseParamsSceneConstrained_splitsSD(:,3),'-','Color',[0 0 0.7])
plot(temporalSupport,pupilData.ellipseParamsSceneConstrained_mean(:,3)+pupilData.ellipseParamsSceneConstrained_splitsSD(:,3),'-','Color',[0 0 0.7])
plot(temporalSupport,pupilData.ellipseParamsAreaSmoothed_mean(:,3),'-r','LineWidth',2)
xlim([0 max(temporalSupport)]);
xlabel('time [mins]');
ylabel('pupil area [pixels in plane]');
hold off

subplot(2,1,2)
plot(temporalSupport,pupilData.meta.smoothPupilArea.pupilAreaMean,'-k');
hold on
plot(temporalSupport,pupilData.meta.smoothPupilArea.pupilAreaMean-pupilData.meta.smoothPupilArea.pupilAreaSD,'-','Color',[0 0 0.7])
plot(temporalSupport,pupilData.meta.smoothPupilArea.pupilAreaMean+pupilData.meta.smoothPupilArea.pupilAreaSD,'-','Color',[0 0 0.7])
plot(temporalSupport,pupilData.meta.smoothPupilArea.pupilAreaMean,'-r','LineWidth',2)
xlim([0 max(temporalSupport)]);
xlabel('time [mins]');
ylabel('pupil area [pixels on eye]');
hold off

% Make a plot of X and Y eye pupil position on the image plane
figure
subplot(2,1,1)
hold on
plot(temporalSupport,pupilData.ellipseParamsUnconstrained_mean(:,1),'-k','LineWidth',2);
hold on
plot(temporalSupport,pupilData.ellipseParamsSceneConstrained_mean(:,1),'-b');
plot(temporalSupport,pupilData.ellipseParamsSceneConstrained_mean(:,1)-pupilData.ellipseParamsSceneConstrained_splitsSD(:,1),'-','Color',[0 0 0.7])
plot(temporalSupport,pupilData.ellipseParamsSceneConstrained_mean(:,1)+pupilData.ellipseParamsSceneConstrained_splitsSD(:,1),'-','Color',[0 0 0.7])
plot(temporalSupport,pupilData.ellipseParamsAreaSmoothed_mean(:,1),'-r','LineWidth',2)
xlim([0 max(temporalSupport)]);
xlabel('time [mins]');
ylabel('X position [pixels]');
hold off

subplot(2,1,2)
hold on
plot(temporalSupport,pupilData.ellipseParamsUnconstrained_mean(:,2),'-k','LineWidth',2);
hold on
plot(temporalSupport,pupilData.ellipseParamsSceneConstrained_mean(:,2),'-b');
plot(temporalSupport,pupilData.ellipseParamsSceneConstrained_mean(:,2)-pupilData.ellipseParamsSceneConstrained_splitsSD(:,2),'-','Color',[0 0 0.7])
plot(temporalSupport,pupilData.ellipseParamsSceneConstrained_mean(:,2)+pupilData.ellipseParamsSceneConstrained_splitsSD(:,2),'-','Color',[0 0 0.7])
plot(temporalSupport,pupilData.ellipseParamsAreaSmoothed_mean(:,2),'-r','LineWidth',2)
xlim([0 max(temporalSupport)]);
xlabel('time [mins]');
ylabel('Y position [pixels]');
hold off
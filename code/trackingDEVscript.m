% tracking demo starting from deinterlaced video

%% clear workspace
clear all
close all
clc

%% Define RUN PARAMS and CONTROL PARAMS
% I think we are at a point where the params struct is being overused, as
% it now includes any possible kind of field one might need throughout the
% entire processing.

% As we move to a more modular code structure, this might lead to errors
% and confusion, as even the simplest function would receive a massive
% param struct as input variable.

% I suggest we keep using the "params strategy" for metadata-kind of
% information (or RUN PARAMS), such as: subject name, session, runName... 

% All tracking parameters (or CONTROL PARAMS) will be fed through an input
% parser into the tracking functions in form of "options" instead. This
% will allow for easier control of each option (manually or via a control
% file), easier default values settings, and it is very much in style with
% matlab's native functions input managment.

%% define RUN PARAMS
% ideally, these are the only necessary standard params to execute the
% whole pupil analysis. Note that this param format is kept almost
% identical also for MRI analysis for code homogeneity.

dropboxDir = '/Volumes/Bay_2_data/giulia/Dropbox-Aguirre-Brainard-Lab/';
params.projectFolder = 'TOME_processing';
params.projectSubfolder = 'session2_spatialStimuli';
params.eyeTrackingDir = 'EyeTracking';

params.subjectName = 'TOME_3020';
params.sessionDate = '050517';
params.runName = 'tfMRI_FLASH_AP_run01';


%% Deinterlace video
% I skip executing this step for now. In principle it would be:

% deinterlaceVideo (params, dropboxDir, 'Mean')


%% %%%%%%%%%%%% TESTING BEGINS HERE %% %%%%%%%%%%%%%%%%%

%% build the input video path

% note that this is the default output format for deinterlaced videos.
inputVideo = fullfile(dropboxDir,params.projectFolder, params.projectSubfolder, ...
        params.subjectName,params.sessionDate,params.eyeTrackingDir, ...
        [params.runName '_60hz.avi']); 

%% prepare the video
disp('Preparing video')

tic
[grayI] = prepareVideo(inputVideo, 'numberOfFrames',500); %just tracking a small portion for testing
toc

% EXEC.TIME = 500 frames take approx 9 sec

%% track the glint
disp('Tracking glint')

tic
glintFile = '/Users/giulia/Desktop/TEST/glintTEST.mat';
[glint, glintTrackingParams] = trackGlint(grayI, glintFile);
toc

% NOTE: displayTracking is not working properly (need to figure out why!).
% Also, as expected if you display the video, the tracking takes longer to
% process.

% IDEA: glint position is very constant, since it is a reflection from a
% fixed source. We might want to exclude samples too far from the mean XY
% position to get rid of false positive traking results.

% EXEC.TIME = 500 frames take approx 60 sec

%% make pupil perimeter video
disp('Making pupil perimeter video')

tic
perimeterVideo = '/Users/giulia/Desktop/TEST/perimeterTEST.avi';
pupilCircleThresh = 0.06; 
pupilEllipseThresh = 0.96;
perimeterParams = extractPupilPerimeter(grayI, perimeterVideo,'pupilCircleThresh', pupilCircleThresh, 'pupilEllipseThresh', pupilEllipseThresh);
toc

% EXEC.TIME = 500 frames take approx 165 sec

%% COMMENTS SO FAR
% 
%  up to this point the routine produces the following output files,
%  necessary for the subsequent steps:
%  1. pupil perimeter video
%  2. glint tracking file (X,Y position frame by frame)
% 
% 
%  The routine also outputs these structs (currently not saved): 
%  1. glintTrackingParams 
%  2. perimeterParams 
%  They include ALL input necessary
%  to replicate the analysis exactly how it was performed the first time
%  around (including the grayI frameseries that originated from prepareVideo).
%  This means that parsing the structs as inputs for the function
%  trackGlint and extractPupilPerimeter respectively will exactly
%  replicate their outputs.
%  The advantage compared to the "params" strategy is again the modularity:
%  only necessary and unambiguous inputs are fed to each step.


%% blink detection
disp('Finding blinks')

tic
% find the blinks
blinkFrames = findBlinks(glintFile);
toc

% show them on the tracked video (this function is for display only)
showBlinks(blinkFrames,grayI)

% note: blinkFrames is an array containing the index of the frames marked as
% blinks.
%% guess pupil cuts
disp('Computing pupil cuts')

tic
framesToCut = guessPupilCuts(perimeterVideo,glintFile,blinkFrames);
toc


%% make control file
controlFileName = '/Users/giulia/Desktop/TEST/controlFileTEST';

makeControlFile(controlFileName, framesToCut, blinkFrames )
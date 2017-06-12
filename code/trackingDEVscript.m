% tracking demo starting from deinterlaced video
clear all
close all
clc

inputVideo = '/Volumes/Bay_2_data/giulia/Dropbox-Aguirre-Brainard-Lab/TOME_processing/session2_spatialStimuli/TOME_3020/050517/EyeTracking/tfMRI_FLASH_AP_run01_60hz.avi';

%% prepare the video
disp('Preparing video')
tic
[grayI] = prepareVideo(inputVideo, 'numberOfFrames',200);
toc

%% track the glint
disp('Tracking glint')
tic
glintFile = '/Users/giulia/Desktop/TEST/glintTEST.mat';
[glint, glintTrackingParams] = trackGlint(grayI, glintFile);
toc
% note: displayTracking is not working properly

% idea: glint position is very constant, since it is a reflection from a
% fixed source. We might want to exclude samples too far from the mean XY
% position to get rid of false positive traking results.


%% threshold pupil video

% DEMO find distances and COP

%% housekeeping
close all
clear all
%clc

%% load pupil data
pupilFileName = '~/Desktop/eyeTrackingDEMO/TOME_processing/session2_spatialStimuli/TOME_3020/050517/EyeTracking/GazeCal01_pupil.mat';

%% find bestFitCoP
bestFitCoP = findBestFitCoP (pupilFileName);
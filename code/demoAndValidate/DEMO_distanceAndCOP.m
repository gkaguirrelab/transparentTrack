% DEMO find distances and COP

%% housekeeping
close all
clear all
%clc

%% load pupil data
pupilFileName = '~/Desktop/eyeTrackingDEMO/TOME_processing/session2_spatialStimuli/TOME_3020/050517/EyeTracking/GazeCal01_pupil.mat';

%% find sceneGeometry
sceneGeometryFileName = '~/Desktop/eyeTrackingDEMO/TOME_processing/session2_spatialStimuli/TOME_3020/050517/EyeTracking/GazeCal01_sceneGeometry.mat';
sceneGeometry = estimateSceneGeometry (pupilFileName,sceneGeometryFileName);

%% given the scene geometry and ellipses centers, find the eccentricity and theta to be used as fitting constraints
% for the demo, imagine that 1 mm = 10 px and the camera was between 10 and
% 14 cm from the subject eye.
distanceRangeInPX = [10*100 10*140];

% load ellises from the pupilfile
load(pupilFileName)
ellipses = pupilData.pPosteriorMeanTransparent;

% initialize eccentricity and theta
eccentricity = nan(length(ellipses),2);
theta = nan(length(ellipses),1);
for ii = 1:length(ellipses)
    if any(isnan(ellipses(ii)))
        continue
    else
        [eccentricity(ii,:), theta(ii)] = constrainEllipseBySceneGeometry (ellipses(ii,:),sceneGeometryFileName,'distanceFromSceneRangePx',distanceRangeInPX);
    end
end
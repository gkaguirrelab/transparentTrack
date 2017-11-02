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

% load ellises from the pupilfile
load(pupilFileName)
ellipses = pupilData.pInitialFitTransparent;
% initialize eccentricity and theta
eccentricity = nan(length(ellipses),1);
theta = nan(length(ellipses),1);
for ii = 1:length(ellipses)
    if any(isnan(ellipses(ii)))
        continue
    else
        [eccentricity(ii,:), theta(ii)] = constrainEllipseBySceneGeometry (ellipses(ii,:),sceneGeometry);
    end
end

figure
plot(ellipses(:,4),eccentricity,'.')
xlabel('transparentEllipses eccentricities')
ylabel('constrained eccentricities')
figure
plot(ellipses(:,5),theta,'.')
xlabel('transparentEllipses tilt angles')
ylabel('constrained tilt angles')

% %% use the derived eccentricity and thetas to constrain the ellipse fit
% % constrained pupil fit here
% 
% %% project the constrained fitted pupil file back in 3d to derive azi and elevation 
% % need to test this
% [reconstructedPupilAzi, reconstructedPupilEle, reconstructedPupilRadius] = pupilProjection_inv(transparentEllipse,centerOfProjection);
% 
% %% do the smoothing in 3d

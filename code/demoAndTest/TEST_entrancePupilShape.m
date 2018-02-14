%% TEST_entrancePupilShape
% Compare our model to results of Mathur 2013
%
% Description:
%   The visual axis of the eye is displaced from the optical (and pupil)
%   axis. Further, the cornea refracts the image of the pupil, causing it
%   to appear magnified and (depending upon viewing angle) shifted. The
%   appearance of the pupil in the image plane is referred to as the
%   "entrance pupil". The shape and size of the entrance pupil will vary as
%   a function of angle with which the eye is viewed. This function has
%   been measured by Spring & Stiles (1948) and Jay (1962), and more
%   recently by Mathur and colleagues:
%
%       Mathur, Ankit, Julia Gehrmann, and David A. Atchison. "Pupil shape
%       as viewed along the horizontal visual field." Journal of vision
%       13.6 (2013): 3-3.
%
%   Given a coordinate system in which viewing angle is calculated with
%   respect to the fixation point of the eye, the entrance pupil will
%   change in appearance as a function of viewing angle differently when
%   viewed from the nasal and temporal side, due to the misalignment of the
%   pupil and visual axes.
%
%   Here, we examine how well our model replicates the empirical work of
%   Mathur 2013.
%

clear all
close all

% Obtain the default sceneGeometry with the following modifications:
% - camera distance to 100 mm
% - center of rotation at the corneal apex
%
% By having the eye rotate around the corneal apex, we match the conditions
% of Mathur et al in which the camera was maintained at a constant
% distance from the corneal apex.
sceneGeometry = createSceneGeometry( ...
    'extrinsicTranslationVector',[0; 0; 100],...
    'eyeLaterality','Right');
sceneGeometry.eye.rotationCenter(1)= 0;

% Assemble the ray tracing functions
rayTraceFuncs = assembleRayTraceFuncs( sceneGeometry );

% Assume a 6 mm true exit pupil diamter, as Mathur 2013 used
% pharmacological dilation for their subjects. The observed entrance pupil
% would have been about 7 mm.
pupilDiam = 6;

% This is Eq 9 from Mathur 2013, which specifies the horizontal to vertical
% ratio of the entrance pupil from different viewing angles relative to
% fixation
mathurEq9 = @(viewingAngleDeg) 0.99.*cosd((viewingAngleDeg+5.3)/1.121);

% This is Eq 11, which specifies the oblique component of pupil ellipticity
mathurEq11 = @(viewingAngleDeg) 0.00072.*viewingAngleDeg-0.0008;

% Mathur 2013 reports results by the visual field angle from which the
% right eye of the subject was observed. A negative angle corresponds to
% viewing the eye from the temporal visual field of the subject.
viewingAngleDeg = -60:1:60;

% Our model rotates the eye. For the right eye, a positive azimuth rotates
% the eye such that the center of the pupil moves to the right of the
% image. This means that a positive azimuth corresponds to the camera being
% positioned in the temporal portion of the visual field. So, we must sign
% reverse the interpretation of our azimuthal values for measurements made
% in the right eye to correspond to the Mathur results. Additionally, we
% need to adjust for kappa: the angle between the pupil and visual axes of
% the eye. The coordinates of our model eye are based around the pupil
% axis. Therfore, we need to calculate a rotation that accounts for the
% Mathur viewing angle and kappa.
azimuthsDeg = (-viewingAngleDeg)-sceneGeometry.eye.kappaAngle(1);
elevationsDeg = zeros(size(viewingAngleDeg))-sceneGeometry.eye.kappaAngle(2);

% Calculate the diameter ratios and thetas
[diamRatios, thetas] = calcPupilDiameterRatio(azimuthsDeg,elevationsDeg,pupilDiam,sceneGeometry,rayTraceFuncs);

% Reverse the thetas to match the Mathur convention, in which a theta of
% zero corresponds to a pupil ellipse with the major axis aligned with the
% horizontal meridian, and positive values of theta are in the
% counter-clockwise direction.
thetas = pi - thetas;

% Calculate the Mathur value C from Equation 6
C = (1-diamRatios).*sin(2.*(thetas-pi/2));

% Plot Figure 10 of Mathur 2013 with our model output.
figure
subplot(1,2,1);
plot(viewingAngleDeg,diamRatios ,'.k');
hold on
plot(viewingAngleDeg,cosd(viewingAngleDeg),'--k');
plot(viewingAngleDeg,mathurEq9(viewingAngleDeg),'-r');
xlim([-90 90]);
ylim([0 1.1]);
xlabel('Viewing angle [deg]')
ylabel('Pupil Diameter Ratio')
title('Mathur 2013 Figure 6, component A')

subplot(1,2,2)
plot(viewingAngleDeg,C ,'.k');
hold on
plot(viewingAngleDeg,mathurEq11(viewingAngleDeg),'-r');
xlim([-90 90]);
ylim([-.2 .2]);
xlabel('Viewing angle [deg]')
ylabel('Oblique component of the pupil ellipticity')
title('Mathur 2013 Figure 6, component C')



%% LOCAL FUNCTION
function [diamRatios, thetas] = calcPupilDiameterRatio(azimuthsDeg,elevationsDeg,pupilDiam,sceneGeometry,rayTraceFuncs)
horizDiam=[];
vertDiam=[];
thetas=[];
for ii = 1:length(azimuthsDeg)
    eyePose=[azimuthsDeg(ii) elevationsDeg(ii) 0 pupilDiam/2];
    % First, perform the forward projection to determine where the center
    % of the pupil is located in the sceneWorld coordinates
    [~, ~, sceneWorldPoints, ~, pointLabels] = pupilProjection_fwd(eyePose, sceneGeometry, rayTraceFuncs, 'fullEyeModelFlag', true);
    % Adjust the sceneGeometry to translate the camera to be centered on
    % center of the pupil. This is not exactly right, as the center of the
    % pupil in the sceneWorld would not exactly correspond to the center of
    % the elliptical entrance pupil seen by the examiner who adjusted the
    % camera, but it is the closest I can get to Atchison's arrangement.
    pupilCenterIdx = find(strcmp(pointLabels,'pupilCenter'));
    adjustedSceneGeometry = sceneGeometry;
    adjustedSceneGeometry.extrinsicTranslationVector(1) = adjustedSceneGeometry.extrinsicTranslationVector(1)+sceneWorldPoints(pupilCenterIdx,1);
    adjustedSceneGeometry.extrinsicTranslationVector(2) = adjustedSceneGeometry.extrinsicTranslationVector(2)+sceneWorldPoints(pupilCenterIdx,2);
    % Now, measure the horizontal and vertical width of the image of the
    % pupil
    [pupilEllipseOnImagePlane, imagePoints] = pupilProjection_fwd(eyePose, adjustedSceneGeometry, rayTraceFuncs, 'nPupilPerimPoints',50);
    horizDiam =[horizDiam max(imagePoints(:,1)')-min(imagePoints(:,1)')];
    vertDiam  =[vertDiam max(imagePoints(:,2)')-min(imagePoints(:,2)')];
    thetas = [thetas, pupilEllipseOnImagePlane(5)];
end
diamRatios=horizDiam./vertDiam;
end



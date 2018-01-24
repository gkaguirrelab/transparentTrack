%% TEST_EntrancePupilMagnification
% Compare pupil magnification in our model to the results of Fedtke 2010
%
% Description:
%   The cornea refracts the image of the pupil, causing it to appear
%   magnified and (depending upon viewing angle) shifted. The appearance of
%   the pupil in the image plane is referred to as the "entrance pupil".
%   Fedtke and colleagues implemented the Navarro schematic model eye in
%   the Zemax ray tracing software package:
%
%       Fedtke, Cathleen, Fabrice Manns, and Arthur Ho. "The entrance pupil
%       of the human eye: a three-dimensional model as a function of
%       viewing angle." Optics express 18.21 (2010): 22364-22376.
%
%   Figure 8 of their paper reports the vertical and horizontal
%   magnification of the entrance pupil, and Figure 9 expresses the ratio
%   of horizontal to vertical entrance pupil size as a function of viewing
%   angle, and compares the output of their model to empirical measurements
%   made by Jay (1962) and Spring & Stiles (1948).
%
%   Here we calculates these values for our model, and compare the results
%   to the prior studies.
%

close all

% Obtain the default sceneGeometry with the following modifications:
% - no lens distortion
% - camera distance to 100 mm
% - center of rotation at the corneal apex

initialCameraPosition = [0; 0; 100];

sceneGeometry = estimateSceneGeometry([],[], ...
    'radialDistortionVector', [0 0], ...
    'extrinsicTranslationVector',initialCameraPosition,...
    'rotationCenterDepth',0);

% Assemble the ray tracing functions
rayTraceFuncs = assembleRayTraceFuncs( sceneGeometry );

% Determine the diameter of a 1.5 mm pupil in the image plane without ray
% tracing effects
pupilDiam = 3;
[~, imagePoints] = pupilProjection_fwd([0 0 0 pupilDiam/2], sceneGeometry, [], 'nPupilPerimPoints',50);
veridicalPupilPixelDiam = max(imagePoints(:,1)')-min(imagePoints(:,1)');

% Rotate the eye between 0 and 80 degrees of azimuth and obtain the
% horizontal and vertical diameter of the pupil image
horizDiam=[];
vertDiam=[];
for azimuthDeg = 0:5:60
    eyeParams=[-azimuthDeg 0 0 pupilDiam/2];
    % We first identify the position of the corneal
    % apex in sceneWorld coordinates following the eye rotation
    [~, ~, sceneWorldPoints, ~, pointLabels] = pupilProjection_fwd(eyeParams, sceneGeometry, rayTraceFuncs, 'fullEyeModelFlag', true);
    cornealApexIdx = find(strcmp(pointLabels,'cornealApex'));
    pupilCenterIdx = find(strcmp(pointLabels,'pupilCenter'));
    % In the Fedtke model, and in the measurements reported by Jay and
    % Spring & Stiles, the camera is positioned such that at translates horizontally to be
    % centered on the center of the pupil.
    tmpSceneGeometry = sceneGeometry;
%    sceneWorldPoints(cornealApexIdx,1:3)
%   tmpSceneGeometry.extrinsicTranslationVector(1) = initialCameraPosition(1) + sceneWorldPoints(pupilCenterIdx,1);
%   tmpSceneGeometry.extrinsicTranslationVector(3) = initialCameraPosition(3) + sceneWorldPoints(pupilCenterIdx,3);
    % Find the image points for this translated camera position.
    [~, imagePoints] = pupilProjection_fwd(eyeParams, tmpSceneGeometry, rayTraceFuncs, 'nPupilPerimPoints',50);
    horizDiam =[horizDiam max(imagePoints(:,1)')-min(imagePoints(:,1)')];
    vertDiam  =[vertDiam max(imagePoints(:,2)')-min(imagePoints(:,2)')];
end

horizMag = horizDiam ./ veridicalPupilPixelDiam;
vertMag = vertDiam ./ veridicalPupilPixelDiam;

% These are equations 1a and 2 of Fedtke 2010
MtanFedtke = @(azimuthDeg, p) (1.133-6.3e-4*p^2)*cosd((-0.8798+4.8e-3*p).*azimuthDeg+3.7e-4.*azimuthDeg.^2);
MsagFedtke = @(azimuthDeg) 4.4e-6.*(azimuthDeg.^2.299)+1.125;

% figure
% subplot(2,1,1);
% plot(0:1:80,MtanFedtke(0:1:80, pupilDiam),'-r');
% hold on
% plot(0:1:80,cosd(0:1:80),'--k');
% plot(0:10:80,horizMag,'xk');
% xlim([0 90]);
% ylim([0.4 1.2]);
% 
% subplot(2,1,2);
% plot(0:1:80,MsagFedtke(0:1:80),'-r');
% hold on
% plot(0:10:80,vertMag,'xk');
% xlim([0 90]);
% ylim([1.1 1.25]);
% 
% figure
% plot(0:1:80,MtanFedtke(0:1:80, pupilDiam)./MsagFedtke(0:1:80),'-r');
% hold on
% plot(0:10:80,horizDiam./vertDiam,'xk');
% xlim([0 90]);
% ylim([0.3 1.1]);

figure
mathurEq = @(azimuthdeg) 0.99.*cosd((azimuthdeg+5.3)/1.121);
plot(-80:1:80,mathurEq(-80:1:80),'-r');
hold on
plot(-80:1:80,cosd(-80:1:80),'--k');
plot(0:-5:-60,horizDiam./vertDiam ,'xk');
plot(0:5:60,horizDiam./vertDiam ,'xk');
xlim([-90 90]);
ylim([0 1.1]);


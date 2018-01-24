%% TEST_EntrancePupilShape
% Compare our model to results of Fedtke 2010 and Mathur 2013
%
% Description:
%   The cornea refracts the image of the pupil, causing it to appear
%   magnified and (depending upon viewing angle) shifted. The appearance of
%   the pupil in the image plane is referred to as the "entrance pupil".
%   The appearance of the entrance pupil as a function of viewing angle has
%   been measured by Spring & Stiles (1948) and Jay (1962), and most
%   recently by Mathur and colleagues:
%
%       Mathur, Ankit, Julia Gehrmann, and David A. Atchison. "Pupil shape
%       as viewed along the horizontal visual field." Journal of vision
%       13.6 (2013): 3-3.
%
%   Fedtke and colleagues created a model of the appearance of the entrance
%   pupil as a function of viewing angle by implementing the Navarro
%   schematic model eye in the Zemax ray tracing software package:
%
%       Fedtke, Cathleen, Fabrice Manns, and Arthur Ho. "The entrance pupil
%       of the human eye: a three-dimensional model as a function of
%       viewing angle." Optics express 18.21 (2010): 22364-22376.
%
%   A key property of the entrance pupil is that it changes in appearance
%   as a function of viewing angle differently when viewed from the nasal
%   and temporal side. This asymmetry is related (at least in part) to the
%   horizontal and vertical decentration of the pupil with respect to the
%   optical axis of the eye.
%
%   Here, we search for values for the center of the pupil relative to the
%   optical axis of the eye that best fits the empirical work of
%   Mathur 2013, and we compare the output of our model to output of the
%   Fedtke model.
%   
%   Note that Mathur 2013 reports results by the visual field angle from
%   which the right eye of the subject was observed. A negative angle
%   corresponds to viewing the eye from the temporal visual field of the
%   subject.
%
%   Our model is rotating the eye. For the right eye, a positive azimuth
%   rotates the eye such that the center of the pupil moves to the right of
%   the image. This means that a positive azimuth corresponds to the camera
%   being positioned in the temporal portion of the visual field. So, we must
%   sign reverse the interpretation of our azimuthal values for measurements
%   made in the right eye to correspond to the Mathur results.
%

close all

% Obtain the default sceneGeometry with the following modifications:
% - no lens distortion
% - camera distance to 100 mm
% - center of rotation at the corneal apex
%
% By having the eye rotate around the corneal apex, we match the conditions
% of these prior papers in which the camera was maintained at a constant
% distance from the corneal apex.
sceneGeometry = estimateSceneGeometry([],[], ...
    'radialDistortionVector', [0 0], ...
    'extrinsicTranslationVector',[0; 0; 100],...
    'rotationCenterDepth',0, ...
    'eyeLaterality','Right');

% Assemble the ray tracing functions
rayTraceFuncs = assembleRayTraceFuncs( sceneGeometry );

% Assume a 6 mm true exit pupil diamter, as Mathur 2013 used
% pharmacological dilation for their subjects. The observed entrance pupil
% would have been about 7 mm.
pupilDiam = 6;

% This is Eq 9 from Mathur 2013, which specifies the horizontal to vertical
% ratio of the entrance pupil from different viewing angles.
mathurEq = @(viewingAngle) 0.99.*cosd((viewingAngle+5.3)/1.121);

% Search over pupil center locations to find the best match to the Marthur
% equation. We first create an objective function that is the difference
% between the horizontal / vertical ratio from our model and Mathur's model
% as a function of viewing angle.
% NOTE: we sign reverse the azimuth here to produce viewing angle.
myObjFunc = @(x) sum((mathurEq(-azimuthsDeg) - calcPupilDiameterRatio(x,azimuthsDeg,pupilDiam,sceneGeometry,rayTraceFuncs)).^2);

% Define some of the parameters of the search.
azimuthsDeg = -60:10:60;
x0=[0.5 -0.25];
ub = [2 2];
lb = [-2 -2];

% Perform the search
[x, fVal] = fmincon(myObjFunc,x0,[],[],[],[],lb,ub);

% Calculate the diameter ratio for the best fitting rotation center values
diamRatio = calcPupilDiameterRatio(x,azimuthsDeg,pupilDiam,sceneGeometry,rayTraceFuncs);

% Plot Figure 10 of Mathur 2013 with our model output
figure
plot(-80:1:80,mathurEq(-80:1:80),'-r');
hold on
plot(-80:1:80,cosd(-80:1:80),'--k');
% NOTE: we sign reverse the azimuth here to produce viewing angle.
plot(-azimuthsDeg,diamRatio ,'xk');
xlim([-90 90]);
ylim([0 1.1]);
xlabel('Viewing angle [deg]')
ylabel('Pupil Diameter Ratio')

% Report the best fitting pupil center
fprintf('For the right eye, the best fitting pupil center was found at:\n');
fprintf('\tp2: %f \n',x(1));
fprintf('\tp3: %f \n',x(2));

%
% horizMag = horizDiam ./ veridicalPupilPixelDiam;
% vertMag = vertDiam ./ veridicalPupilPixelDiam;
%
% % These are equations 1a and 2 of Fedtke 2010
% MtanFedtke = @(azimuthDeg, p) (1.133-6.3e-4*p^2)*cosd((-0.8798+4.8e-3*p).*azimuthDeg+3.7e-4.*azimuthDeg.^2);
% MsagFedtke = @(azimuthDeg) 4.4e-6.*(azimuthDeg.^2.299)+1.125;
%
% figure
% subplot(2,1,1);
% plot(0:1:80,MtanFedtke(0:1:80, pupilDiam),'-r');
% hold on
% plot(0:1:80,cosd(0:1:80),'--k');
% plot(0:5:60,horizMag,'xk');
% xlim([0 90]);
% ylim([0.4 1.2]);
%
% subplot(2,1,2);
% plot(0:1:80,MsagFedtke(0:1:80),'-r');
% hold on
% plot(0:5:60,vertMag,'xk');
% xlim([0 90]);
% ylim([1.1 1.25]);
%
% figure
% plot(0:1:80,MtanFedtke(0:1:80, pupilDiam)./MsagFedtke(0:1:80),'-r');
% hold on
% plot(0:5:60,horizDiam./vertDiam,'xk');
% xlim([0 90]);
% ylim([0.3 1.1]);
%


%% LOCAL FUNCTION
function diamRatio = calcPupilDiameterRatio(x,azimuthsDeg,pupilDiam,sceneGeometry,rayTraceFuncs)
horizDiam=[];
vertDiam=[];
% Update the sceneGeometry with the past center of rotation value
sceneGeometry.eye.pupilCenter = [-3.7,x(1),x(2)];
for ii = 1:length(azimuthsDeg)
    eyeParams=[azimuthsDeg(ii) 0 0 pupilDiam/2];
    % First, perform the forward projection to determine where the center
    % of the pupil is located in the sceneWorld coordinates
    [~, ~, sceneWorldPoints, ~, pointLabels] = pupilProjection_fwd(eyeParams, sceneGeometry, rayTraceFuncs, 'fullEyeModelFlag', true);
    % Adjust the sceneGeometry to translate the camera to be centered on
    % center of the pupil. This is not exactly right, as the center of the
    % pupil in the sceneWorld would not exactly correspond to the center of
    % the elliptical entrance pupil seen by the examiner who adjusted the
    % camera, but it is the closest I can get to Atchison's arrangement.
    pupilCenterIdx = find(strcmp(pointLabels,'pupilCenter'));
    sceneGeometry.extrinsicTranslationVector(1) = sceneGeometry.extrinsicTranslationVector(1)+sceneWorldPoints(pupilCenterIdx,1);
    % Now, measure the horizontal and vertical width of the image of the
    % pupil
    [~, imagePoints] = pupilProjection_fwd(eyeParams, sceneGeometry, rayTraceFuncs, 'nPupilPerimPoints',50);
    horizDiam =[horizDiam max(imagePoints(:,1)')-min(imagePoints(:,1)')];
    vertDiam  =[vertDiam max(imagePoints(:,2)')-min(imagePoints(:,2)')];
end
diamRatio=horizDiam./vertDiam;
end



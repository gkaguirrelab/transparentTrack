%% TEST_entrancePupilShape
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
%   A key property of the entrance pupil is that it changes in appearance
%   as a function of viewing angle differently when viewed from the nasal
%   and temporal side. This asymmetry is related (at least in part) to the
%   horizontal and vertical decentration of the pupil with respect to the
%   optical axis of the eye.
%
%   Here, we search for values for the center of the pupil relative to the
%   optical axis of the eye that best fits the empirical work of
%   Mathur 2013.
%   
%   Note that Mathur 2013 reports results by the visual field angle from
%   which the right eye of the subject was observed. A negative angle
%   corresponds to viewing the eye from the temporal visual field of the
%   subject.
%
%   Our model rotates the eye. For the right eye, a positive azimuth
%   rotates the eye such that the center of the pupil moves to the right of
%   the image. This means that a positive azimuth corresponds to the camera
%   being positioned in the temporal portion of the visual field. So, we
%   must sign reverse the interpretation of our azimuthal values for
%   measurements made in the right eye to correspond to the Mathur results.
%

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
% ratio of the entrance pupil from different viewing angles.
mathurEq = @(viewingAngle) 0.99.*cosd((viewingAngle+5.3)/1.121);

% Search over pupil center locations to find the best match to the Marthur
% equation.

% Define some of the parameters of the search.
azimuthsDeg = -60:10:60;
x0=sceneGeometry.eye.pupilCenter(2:3);
% As the solution is symmetric for p3 values around zero, we make the lower
% bound on the p3 value zero to place the resulting pupil center downward
% from the corneal apex.
lb = [-1 0];
ub = [1 1];

% Create an objective function that is the difference
% between the horizontal / vertical ratio from our model and Mathur's model
% as a function of viewing angle.
% NOTE: we sign reverse the azimuth here to produce viewing angle.
myObjFunc = @(x) sum((mathurEq(-azimuthsDeg) - calcPupilDiameterRatio(x,azimuthsDeg,pupilDiam,sceneGeometry,rayTraceFuncs)).^2);

% Perform the search
[x, fVal] = fmincon(myObjFunc,x0,[],[],[],[],lb,ub);

% Calculate the diameter ratio for the best fitting rotation center values
diamRatio = calcPupilDiameterRatio(x,azimuthsDeg,pupilDiam,sceneGeometry,rayTraceFuncs);

% Plot Figure 10 of Mathur 2013 with our model output
figure
plot(-azimuthsDeg,mathurEq(-azimuthsDeg),'-r');
hold on
plot(-azimuthsDeg,cosd(-azimuthsDeg),'--k');
% NOTE: plot the diamRatio against viewingAngle
plot(azimuthsDeg,diamRatio ,'xk');
xlim([-90 90]);
ylim([0 1.1]);
xlabel('Viewing angle [deg]')
ylabel('Pupil Diameter Ratio')

% Report the best fitting pupil center
fprintf('For the right eye, the best fitting pupil center was found at:\n');
fprintf('\tp2: %f \n',x(1));
fprintf('\tp3: %f \n',x(2));




%% LOCAL FUNCTION
function diamRatio = calcPupilDiameterRatio(x,azimuthsDeg,pupilDiam,sceneGeometry,rayTraceFuncs)
horizDiam=[];
vertDiam=[];
% Update the sceneGeometry with the past center of rotation value
sceneGeometry.eye.pupilCenter(2:3) = x;
for ii = 1:length(azimuthsDeg)
    eyePose=[azimuthsDeg(ii) 0 0 pupilDiam/2];
    % First, perform the forward projection to determine where the center
    % of the pupil is located in the sceneWorld coordinates
    [~, ~, sceneWorldPoints, ~, pointLabels] = pupilProjection_fwd(eyePose, sceneGeometry, rayTraceFuncs, 'fullEyeModelFlag', true);
    % Adjust the sceneGeometry to translate the camera to be centered on
    % center of the pupil. This is not exactly right, as the center of the
    % pupil in the sceneWorld would not exactly correspond to the center of
    % the elliptical entrance pupil seen by the examiner who adjusted the
    % camera, but it is the closest I can get to Atchison's arrangement.
    pupilCenterIdx = find(strcmp(pointLabels,'pupilCenter'));
    sceneGeometry.extrinsicTranslationVector(1) = sceneGeometry.extrinsicTranslationVector(1)+sceneWorldPoints(pupilCenterIdx,1);
    % Now, measure the horizontal and vertical width of the image of the
    % pupil
    [~, imagePoints] = pupilProjection_fwd(eyePose, sceneGeometry, rayTraceFuncs, 'nPupilPerimPoints',50);
    horizDiam =[horizDiam max(imagePoints(:,1)')-min(imagePoints(:,1)')];
    vertDiam  =[vertDiam max(imagePoints(:,2)')-min(imagePoints(:,2)')];
end
diamRatio=horizDiam./vertDiam;
end



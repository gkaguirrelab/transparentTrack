%% DEMO_eyePoseParams
% Demonstrate the interpretation of the eyePose values
%
% Description:
%   Our routines express the pose of the eye in a 1x4 vector of eyePoses
%   of the form [azimuth, elevation, torsion, pupilRadius]
%
% Here we demonstrate the interpretation of these parameters
%

clear all
close all
clc

fprintf(['The pose of the eye is described by the parameters:\n\n' ...
    '\t[azimuth, elevation, torsion, pupilRadius]\n\n' ...
    'The three rotation variables are in units of degrees, and are in the\n' ...
    'head-fixed, extrinsic coordinate space. This means that the parameters\n' ...
    'are unlike "Fick" coordinates, which are with reference to an intrinsic\n' ...
    '(i.e., rotating) coordinate frame.\n\n']);

% Obtain the sceneGeometry and ray tracing functions
sceneGeometry = createSceneGeometry('eyeLaterality','Right');
rayTraceFuncs = assembleRayTraceFuncs( sceneGeometry );

% Define some variables for plotting the model eye
eyePartLabels = {'rotationCenter', 'posteriorChamber' 'irisPerimeter' 'pupilPerimeter' 'anteriorChamber' 'cornealApex' 'pupilCenter'};
plotColors = {'+r' '.w' '.b' '*g' '.y' '*y' '+g'};
blankFrame = zeros(480,640)+0.5;

%% Present Figure 1
figure(1)
eyePoses=[-20 20 0 3; 0 20 0 3; 20 20 0 3; -20 0 0 3; 0 0 0 3; 20 0 0 3; -20 -20 0 3; 0 -20 0 3; 20 -20 0 3 ];

for pose = 1:size(eyePoses,1)
    % Perform the projection and request the full eye model
    [~, imagePoints, ~, ~, pointLabels] = pupilProjection_fwd(eyePoses(pose,:),sceneGeometry,rayTraceFuncs,'fullEyeModelFlag',true);
    % plot
    subplot(3,3,pose);
    imshow(blankFrame, 'Border', 'tight');
    hold on
    axis off
    axis equal
    xlim([0 640]);
    ylim([0 480]);
    % Plot each anatomical component
    for pp = 1:length(eyePartLabels)-1
        idx = strcmp(pointLabels,eyePartLabels{pp});
        plot(imagePoints(idx,1), imagePoints(idx,2), plotColors{pp})
    end
    title(num2str(eyePoses(pose,:)));
end
drawnow
fprintf(['Figure 1 shows the pose of the eye across positive and negative values\n' ...
    'of azimuth and elevation. Torsion is fixed at zero. For saccadic\n' ...
    'eye movements under head-fixed conditions, Listing`s Law holds that \n' ...
    'rotations of the eye keep torsion constant.\n\n']);

%% Present Figure 2
figure(2)
eyeSides = {'right','left'};
for laterality = 1:2
    % prepare the model eye for this laterality
    sceneGeometry = createSceneGeometry('eyeLaterality',eyeSides{laterality});
    rayTraceFuncs = assembleRayTraceFuncs( sceneGeometry );
    [pupilEllipseOnImagePlane, imagePoints, ~, ~, pointLabels] = pupilProjection_fwd([0 0 0 3],sceneGeometry,rayTraceFuncs,'fullEyeModelFlag',true);
    
    % setup the figure
    subplot(1,2,laterality);
    imshow(blankFrame, 'Border', 'tight');
    hold on
    axis off
    axis equal
    xlim([0 640]);
    ylim([0 480]);
    text(320,50,[eyeSides{laterality} ' eye'],'HorizontalAlignment','center');
    
    % Plot each anatomical component
    pFitImplicit = ellipse_ex2im(ellipse_transparent2ex(pupilEllipseOnImagePlane));
    fh=@(x,y) pFitImplicit(1).*x.^2 +pFitImplicit(2).*x.*y +pFitImplicit(3).*y.^2 +pFitImplicit(4).*x +pFitImplicit(5).*y +pFitImplicit(6);
    fimplicit(fh,[1, 640, 1, 480],'Color', 'g','LineWidth',1);
    axis off;
    idx = strcmp(pointLabels,eyePartLabels{4});
    plot(imagePoints(idx,1), imagePoints(idx,2), plotColors{4})
    idx = strcmp(pointLabels,eyePartLabels{6});
    plot(imagePoints(idx,1), imagePoints(idx,2), plotColors{6})
    idx = strcmp(pointLabels,eyePartLabels{7});
    plot(imagePoints(idx,1), imagePoints(idx,2), plotColors{7})
    hold off
end
drawnow
fprintf(['Figure 2 shows just the perimeter of the pupil and the corneal apex\n' ...
    'for eyePoses [0 0 0 3] for the right and left eye. Here, the axis of\n' ...
    'the camera is aligned with the optical axis of the model eye. \n' ...
    'Note that the center of the pupil is displaced downwards and nasally\n' ...
    'with respect to the optical axis of each eye. This physiologic property\n' ...
    'causes the entrance pupil to have a slightly different appearance\n' ...
    'when viewed from the nasal or temporal visual field (Atchison, 2013).\n\n']);

%% Present Figure 3
figure(3)
sceneGeometry = createSceneGeometry();
rayTraceFuncs = assembleRayTraceFuncs( sceneGeometry );

% Plot the values for eyePost [0 0 0 x]
subplot(1,2,1);
entrancePupilRadiusRayTrace = [];
entrancePupilRadiusNoRayTrace = [];
for radius = 0.5:0.5:4
    pupilEllipseOnImagePlane = pupilProjection_fwd([0 0 0 radius],sceneGeometry,rayTraceFuncs,'fullEyeModelFlag',true);
    entrancePupilRadiusRayTrace = [entrancePupilRadiusRayTrace sqrt(pupilEllipseOnImagePlane(3)/pi)];
    pupilEllipseOnImagePlane = pupilProjection_fwd([0 0 0 radius],sceneGeometry,[],'fullEyeModelFlag',true);
    entrancePupilRadiusNoRayTrace = [entrancePupilRadiusNoRayTrace sqrt(pupilEllipseOnImagePlane(3)/pi)];
end
plot(0.5:0.5:4,entrancePupilRadiusNoRayTrace,'*k');
hold on
plot(0.5:0.5:4,entrancePupilRadiusRayTrace,'*r');
lsline()
hold off
xlim([0 5]);
ylim([0 125]);
xlabel('pupil radius modeled [mm]');
ylabel('entrance pupil radius in image [pixels]');
legend({'Without corneal refraction','With corneal refraction'},'Location','southeast');
title('eyePose [0 0 0 x]');

% Plot the values for eyePost [15 15 0 x]
subplot(1,2,2);
entrancePupilRadiusRayTrace = [];
entrancePupilRadiusNoRayTrace = [];
for radius = 0.5:0.5:4
    pupilEllipseOnImagePlane = pupilProjection_fwd([15 15 0 radius],sceneGeometry,rayTraceFuncs,'fullEyeModelFlag',true);
    entrancePupilRadiusRayTrace = [entrancePupilRadiusRayTrace sqrt(pupilEllipseOnImagePlane(3)/pi)];
    pupilEllipseOnImagePlane = pupilProjection_fwd([15 15 0 radius],sceneGeometry,[],'fullEyeModelFlag',true);
    entrancePupilRadiusNoRayTrace = [entrancePupilRadiusNoRayTrace sqrt(pupilEllipseOnImagePlane(3)/pi)];
end
plot(0.5:0.5:4,entrancePupilRadiusNoRayTrace,'*k');
hold on
plot(0.5:0.5:4,entrancePupilRadiusRayTrace,'*r');
lsline()
hold off
xlim([0 5]);
ylim([0 125]);
xlabel('pupil radius modeled [mm]');
ylabel('entrance pupil radius in image [pixels]');
legend({'Without corneal refraction','With corneal refraction'},'Location','southeast');
title('eyePose [30 25 0 x]');
drawnow

fprintf(['Figure 3 shows the radius of the entrance pupil in the image\n' ...
    'as a function of the modeled radius of the physical pupil for an\n' ...
    'eye in two different poses. The black points are the values for \n' ...
    'a model in which the refractive effects of the cornea are not modeled \n' ...
    'while the red points do include this component. Corneal refraction \n' ...
    'magnifies the image of the pupil and, depending upon eye rotation, \n' ...
    'will also alter the elliptical shape of the pupil on the image plane.\n']);

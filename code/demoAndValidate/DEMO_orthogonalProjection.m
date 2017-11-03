% DEMO_orthogonalProjection
% this demo explores the problem of orthogonal projection

% current assumptions:
% 1. the pupil is a circle of constant radius living on a plane that can
% rotate freely around an axis of length R and centered in the origin.
% 2. the scene is a 2D plane at a known position from the center of
% rotation.
% 3. the center of the scene is the center of projection. When the
% projection axis is coincident with the pupil plane normal axis, the
% projected pupil on the scene plane will have a circular point.
% 4. units for angles are degrees, linear units are uniform and arbitrary.

%% housekeeping
clear all
close all
clc

%% construct initial geometry

% to obtain assumption 1, we intersect a sphere with a radius slightly
% bigger than the pupil plane rotation axis with a plane orthogonal with
% the surface of the sphere at a distance R from the origin. The
% intersection will result in a circle of constant radius on the pupil
% plane, no matter the orientation of the plane itself. NOTE that the
% intersecting sphere by which we define the pupil plane at this point is
% not meant to represent the eyeball. It is merely a "trick" to make sure
% that the pupil is always a circle of the same size rotating around a
% fixed arm in 3D.

% define the eyeball sphere radius
eyeballR =  100; 
eyeballCenter = [10 0 0];
eyeball = [eyeballCenter eyeballR];

% "Depth" of the pupil plane, with respect to the sphere radius. The
% deeper, the bigger the pupil.
planeDepth = 5;
intersectingSphere = [eyeballCenter eyeballR-planeDepth];


% "Depth" of the intersecting plane, with respect to the sphere radius. The
% deeper, the bigger the pupil.
planeDepth = 5;

% rotation arm length
rotationArmLength = eyeballR - planeDepth;

% create scene plane
sceneDistance = 100; % orthogonal distance from rotation arm
xMax= 320/2; % max x size of scene (for plotting purposes)
yMax = 240/2; % max y size of scene
scenePlane = createPlane([0 0 rotationArmLength+sceneDistance],[0 0 rotationArmLength+sceneDistance]);

% center of projection
centerOfProjection3D = projPointOnPlane(eyeballCenter, scenePlane);


% derive optical axis
% [centerOfProjectionX,centerOfProjectionY,centerOfProjectionZ] = sph2cart(azi0,ele0,rotationArmLength);
opticalAxis = createLine3d(eyeballCenter,centerOfProjection3D);

%% rotate the optical axis along X and Y and draw the pupil
% define rotations in deg
pupilAzi = 20; % in degrees
pupilEle = 0; % in degrees

% first rotate along Y (azi)
rotationX = createRotationOy(eyeballCenter, deg2rad(pupilAzi));
pupilAxis = transformLine3d(opticalAxis, rotationX);

% then along X (elevation)
rotationY = createRotationOx(eyeballCenter, deg2rad(-pupilEle)); % - to preserve convention
pupilAxis = transformLine3d(pupilAxis, rotationY);

% intersect sphere and pupil axis
pupilCenter3D = intersectLineSphere(pupilAxis, intersectingSphere)

% create the pupil plane
pupilPlane = createPlane(pupilCenter3D(2,:),pupilCenter3D(2,:)-eyeballCenter);

% intersect pupil plane and sphere to obtain the pupil circle
pupilInEye = intersectPlaneSphere(pupilPlane, eyeball);
pupilRadius = pupilInEye(4); % saving it explicitely for reference

%  get the 3d cartesian coordinates of the pupil perimeter points
pupilPoints3d = getCirclePoints3d(pupilInEye);

% do orthogonal projeciton on the scene plane
pupilPoints2d = projPointOnPlane(pupilPoints3d,scenePlane);

%% do ellipse fitting in scene plane
% uncostrained ellipse fit
unconstrainedEllipseImplicit = ellipsefit_direct(pupilPoints2d(:,1), pupilPoints2d(:,2));
unconstrainedEllipseTransparent = ellipse_ex2transparent(ellipse_im2ex(unconstrainedEllipseImplicit))
unconstrainedEllipseExplicit = ellipse_im2ex(unconstrainedEllipseImplicit);

%% Test projections functions
 % define center of projection (ortogonal approx => no perspective
 % displacement from 3D to 2D)
centerOfProjection = [centerOfProjection3D(1) centerOfProjection3D(2)];

% inverse projection
[reconstructedPupilAzi, reconstructedPupilEle, reconstructedPupilRadius] = pupilProjection_inv(unconstrainedEllipseTransparent,centerOfProjection);

% forward projection
reconstructedTransparentEllipse = pupilProjection_fwd(pupilAzi, pupilEle, pupilCenter3D(2,:), 'pupilRadius', pupilRadius);


%% test constraints on eccentricity and theta
% find candidates eyeball
candidatesEB = findCandidatesEyeball(reconstructedTransparentEllipse,eyeballR-planeDepth);

% since we have just a single ellipse, we just pick the first candidate
% eyeball 

sceneGeometry.eyeballCenter.X = candidatesEB(1,1);
sceneGeometry.eyeballCenter.Y = candidatesEB(1,2);
sceneGeometry.eyeballCenter.Z = candidatesEB(1,3);% meant as the distance from the scene plane
sceneGeometry.eyeballRadius = eyeballR-planeDepth;


% reconstruct theta and eccentricity
[eccentricity, theta] = constrainEllipseBySceneGeometry (reconstructedTransparentEllipse,sceneGeometry);

figure
plot(reconstructedTransparentEllipse(4),eccentricity,'.')
xlabel('transparentEllipse eccentricities')
ylabel('constrained eccentricities')
figure
plot(reconstructedTransparentEllipse(5),theta,'.')
xlabel('transparentEllipses tilt angles')
ylabel('constraine tilt angles')

%%  plots

% plot 3d
figure
grid on

% add the scene
patch([xMax -xMax -xMax xMax], [yMax yMax -yMax -yMax], [rotationArmLength+sceneDistance rotationArmLength+sceneDistance rotationArmLength+sceneDistance rotationArmLength+sceneDistance], 'yellow','FaceAlpha',0.5)
hold on

% add the eyeball
drawSphere(eyeball,'FaceAlpha',0.5)
hold on

% add optical axis
drawLine3d(opticalAxis,'k')
hold on

% add the pupil on the eye
drawCircle3d(pupilInEye)
hold on

% add pupil axis
drawLine3d(pupilAxis,'b')
hold on

% add the pupil to the scene in the figure
plot3(pupilPoints2d(:,1),pupilPoints2d(:,2),rotationArmLength+sceneDistance *ones(length(pupilPoints2d)),'-','color','red')

% add labels
xlabel('X')
ylabel('Y')
zlabel('Z')
title( 'Forward orthogonal projection (from eyeball to scene)')
axis equal

% plot 2d scene
figure

% plot original points
plot(pupilPoints2d(:,1),pupilPoints2d(:,2))
hold on

% plot transparent ellipse
[eX, eY, ~, ~] = ellipse(500,unconstrainedEllipseExplicit(1),unconstrainedEllipseExplicit(2),unconstrainedEllipseExplicit(3),unconstrainedEllipseExplicit(4),unconstrainedEllipseExplicit(5));
plot(eX, eY,'.')
hold on
% plot reconstructed ellipse
reconstructedEllipseExplicit = ellipse_transparent2ex(reconstructedTransparentEllipse);
[rX, rY, ~, ~] = ellipse(500,reconstructedEllipseExplicit(1),reconstructedEllipseExplicit(2),reconstructedEllipseExplicit(3),reconstructedEllipseExplicit(4),reconstructedEllipseExplicit(5));
plot(rX, rY, '.')
hold off
legend('Projection', 'Ellipse Fit', 'Reconstruction')
xlim([-xMax xMax])
ylim([-yMax yMax])
title('View of the scene')
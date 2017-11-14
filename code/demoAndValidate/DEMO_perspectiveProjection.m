% DEMO_perspectiveProjection
% this demo explores the problem of orthogonal projection

% current assumptions:
% 1. the pupil is a circle of constant radius living on a plane that can
% rotate freely around an axis of length R and centered in the origin.
% 2. the scene is a 2D plane at a known position from the center of
% rotation.
% 3. units for angles are degrees, linear units are uniform and arbitrary.

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
% not meant to represent the eye. It is merely a "trick" to make sure
% that the pupil is always a circle of the same size rotating around a
% fixed arm in 3D.

% define the eye sphere radius
eyeR =  250; 
eyeCenter = [0 0 0];
eye = [eyeCenter eyeR];


% "Depth" of the pupil plane, with respect to the sphere radius. The
% deeper, the bigger the pupil.
planeDepth = 5;
intersectingSphere = [eyeCenter eyeR-planeDepth];


% rotation arm length
rotationArmLength = eyeR - planeDepth;

% create scene plane
sceneDistance = 1200; % orthogonal distance from rotation arm
xMax= 320; % max x size of scene (for plotting purposes)
yMax = 240; % max y size of scene
scenePlane = createPlane([0 0 rotationArmLength+sceneDistance],[0 0 rotationArmLength+sceneDistance]);

% center of projection
centerOfProjection3D = projPointOnPlane(eyeCenter, scenePlane);


% derive optical axis
% [centerOfProjectionX,centerOfProjectionY,centerOfProjectionZ] = sph2cart(azi0,ele0,rotationArmLength);
opticalAxis = createLine3d(eyeCenter,centerOfProjection3D);


%% rotate the optical axis along X and Y and draw the pupil
% define rotations in deg
pupilAzi = 20; % in degrees
pupilEle = 0; % in degrees

% first rotate along Y (azi)
rotationX = createRotationOy(eyeCenter, deg2rad(pupilAzi));
pupilAxis = transformLine3d(opticalAxis, rotationX);

% then along X (elevation)
rotationY = createRotationOx(eyeCenter, deg2rad(-pupilEle)); % - to preserve convention
pupilAxis = transformLine3d(pupilAxis, rotationY);

% intersect sphere and pupil axis
pupilCenter3D = intersectLineSphere(pupilAxis, intersectingSphere);

% create the pupil plane
pupilPlane = createPlane(pupilCenter3D(2,:),pupilCenter3D(2,:)-eyeCenter);

% intersect pupil plane and sphere to obtain the pupil circle
pupilInEye = intersectPlaneSphere(pupilPlane, eye);
pupilRadius = pupilInEye(4); % saving it explicitely for reference

%  get the 3d cartesian coordinates of the pupil perimeter points
pupilPoints3d = getCirclePoints3d(pupilInEye);


%% reconstruct scene geometry.




%%  plots

% plot 3d
figure
grid on

% add the scene
patch([xMax -xMax -xMax xMax], [yMax yMax -yMax -yMax], [rotationArmLength+sceneDistance rotationArmLength+sceneDistance rotationArmLength+sceneDistance rotationArmLength+sceneDistance], 'yellow','FaceAlpha',0.5)
hold on

% add the eyeball
drawSphere(eye,'FaceAlpha',0.5)
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


% add labels
xlabel('X')
ylabel('Y')
zlabel('Z')

axis equal

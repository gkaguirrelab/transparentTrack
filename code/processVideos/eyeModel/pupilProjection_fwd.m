function [pupilEllipseOnImagePlane, virtualEyeWorldPoints, imagePoints, pointLabels] = pupilProjection_fwd(eyeParams, sceneGeometry, corneaRayTraceFunc, fullEyeModel)
% Project the pupil circle to an ellipse on the image plane
%
% Description:
%	Given the sceneGeometry, this routine simulates a circular pupil on a
%   spherical eye and then measures the parameters of the ellipse (in
%   transparent format) of the projection of the pupil to the image plane.
%
% Notes:
%   Rotations - Eye rotations are given as azimuth and elevations in
%   degrees. These values correspond to degrees of rotation of the eye
%   relative to a head-fixed (extrinsic) coordinate frame. Note that this
%   is different from an eye-fixed (intrinsic) coordinate frame (such as
%   the Fick coordinate sysem). Azimuth, Elevation of [0,0] corresponds
%   to the position of the eye when a line that connects the center of
%   rotation of the eye with the center of the pupil is normal to the image
%   plane. Positive rotations correspond to rightward, upward, translation
%   of the pupil center in the image.
%
%   Units - Eye rotations are in units of degrees. However, the units of
%   theta in the transparent ellipse parameters are radians. This is in
%   part to help us keep the two units separate conceptually.
%
% Inputs:
%   eyeParams             - A 1x3 vector provides values for [eyeAzimuth,
%                           eyeElevation, pupilRadius]. Azimuth and
%                           elevation are in units of head-centered
%                           (extrinsic) degrees, and pupil radius is in mm.
%   sceneGeometry         - A structure that contains the fields:
%                             - A 1x2 vector of radial distortion params
%                             - intrinsicCameraMatrix: a 3x3 matrix in
%                               arbitrary units (typically pixels)
%                             - extrinsicTranslationVector: a 3x1 vector
%                               in units of millimeters that relates center
%                               of rotation of the eye to the optical axis
%                               of the camera
%                             - extrinsicRotationMatrix: a 3x3 matrix in
%                               units of millimeters
%                             - eye: a sub-structure that contains fields
%                               that define the anatomical properties of
%                               the eye. The values are largely those
%                               returned by modelEyeParameters(), with the
%                               rotationCenter field adjusted by the
%                               sceneGeometry search
%   fullEyeModel          - Logical. Determines if the full posterior and
%                           anterior chamber eye model will be created.
%
% Outputs:
%   pupilEllipseOnImagePlane - A 1x5 vector that contains the parameters of
%                           pupil ellipse on the image plane cast in
%                           transparent form
%   eyeWorldPoints        - An nx3 matrix that gives the coordinates of the
%                           points of the eye model in the eyeWorld
%                           coordinate frame. If fullEyeModel is set to
%                           false, then n=5 for just the pupil perimeter.
%                           If fullEyeModel is true, then n ~500.
%   imagePoints           - An nx2 matrix that specifies the x, y location
%                           on the image plane for each of the eyeWorld
%                           points.
%   pointsLabels          - An nx1 cell array that identifies each of the
%                           points, from the set {'pupilCenter',
%                           'irisCenter', 'rotationCenter',
%                           'posteriorChamber', 'irisPerimeter',
%                           'pupilPerimeter', 'anteriorChamber'}.
%

%% Check the input
if nargin==0
    % No inputs were provided
    error('Provide the parameters for the pose of the eye [eyeAzimuth, eyeElevation, pupilRadius].');
end
if nargin==1
    % No sceneGeometry was provided. Use the default settings
    sceneGeometry = estimateSceneGeometry('','');
end
if nargin<=2
    % No ray trace function was provided. Set to empty
    corneaRayTraceFunc = [];
    corneaRayTraceFunc = createCorneaRayTraceFunction( sceneGeometry );
end
if nargin<=3
    % No valye was passed for the fullEyeModel flag. Set to false.
    fullEyeModel = false;
end


%% Prepare variables
% Separate the eyeParams into individual variables
eyeAzimuth = eyeParams(1);
eyeElevation = eyeParams(2);
pupilRadius = eyeParams(3);


%% Define an eye in eyeWorld coordinates
% This coordinate frame is in mm units and has the dimensions (p1,p2,p3).
% The diagram is of a cartoon pupil, being viewed directly from the front.
%
% Coordinate [0,0,0] corresponds to the apex (front surface) of the cornea.
% The first dimension is depth, and has a negative value towards the
% back of the eye.
%
%                 |
%     ^         __|__
%  +  |        /     \
% p2  -  -----(   +   )-----
%  -  |        \_____/
%     v           |
%                 |
%
%           - <--p3--> +
%


% Five points are defined around the pupil circle, which uniquely
% constrains the ellipse in the image plane.
nPupilPerimPoints = 5;
perimeterPointAngles = 0:2*pi/nPupilPerimPoints:2*pi-(2*pi/nPupilPerimPoints);
eyeWorldPoints(1:nPupilPerimPoints,3) = ...
    sin(perimeterPointAngles)*pupilRadius + sceneGeometry.eye.pupilCenter(3);
eyeWorldPoints(1:nPupilPerimPoints,2) = ...
    cos(perimeterPointAngles)*pupilRadius + sceneGeometry.eye.pupilCenter(2);
eyeWorldPoints(1:nPupilPerimPoints,1) = ...
    0 + sceneGeometry.eye.pupilCenter(1);

% Create labels for the pupilPerimeter points
tmpLabels = cell(nPupilPerimPoints, 1);
tmpLabels(:) = {'pupilPerimeter'};
pointLabels = tmpLabels;

% If the fullEyeModel flag is set, then we will create a set of points that
% define an anatomical model of the posterior and anterior chambers of the
% eye.
if fullEyeModel
    
    % Add points for the center of the pupil, iris, and rotation
    eyeWorldPoints = [eyeWorldPoints; sceneGeometry.eye.pupilCenter];
    pointLabels = [pointLabels; 'pupilCenter'];
    eyeWorldPoints = [eyeWorldPoints; sceneGeometry.eye.irisCenter];
    pointLabels = [pointLabels; 'irisCenter'];
    eyeWorldPoints = [eyeWorldPoints; sceneGeometry.eye.rotationCenter];
    pointLabels = [pointLabels; 'rotationCenter'];
    
    % Create the posterior chamber ellipsoid. We switch dimensions here so
    % that the ellipsoid points have their poles at corneal apex and
    % posterior apex of the eye
    [p3tmp, p2tmp, p1tmp] = ellipsoid( ...
        sceneGeometry.eye.posteriorChamberCenter(3), ...
        sceneGeometry.eye.posteriorChamberCenter(2), ...
        sceneGeometry.eye.posteriorChamberCenter(1), ...
        sceneGeometry.eye.posteriorChamberRadii(3), ...
        sceneGeometry.eye.posteriorChamberRadii(2), ...
        sceneGeometry.eye.posteriorChamberRadii(1), ...
        30);
    % Convert the surface matrices to a vector of points and switch the
    % axes back
    ansTmp = surf2patch(p1tmp, p2tmp, p3tmp);
    posteriorChamberPoints=ansTmp.vertices;
    
    % Retain those points that are anterior to the center of the posterior
    % chamber and are posterior to the iris plane
    retainIdx = logical(...
        (posteriorChamberPoints(:,1) > sceneGeometry.eye.posteriorChamberCenter(1)) .* ...
        (posteriorChamberPoints(:,1) < sceneGeometry.eye.irisCenter(1)) ...
        );
    if all(~retainIdx)
        error('The iris center is behind the center of the posterior chamber');
    end
    posteriorChamberPoints = posteriorChamberPoints(retainIdx,:);
    
    % Add the points and labels
    eyeWorldPoints = [eyeWorldPoints; posteriorChamberPoints];
    tmpLabels = cell(size(posteriorChamberPoints,1), 1);
    tmpLabels(:) = {'posteriorChamber'};
    pointLabels = [pointLabels; tmpLabels];
    
    % Define 360 points around the perimeter of the iris
    nIrisPerimPoints = 360;
    perimeterPointAngles = 0:2*pi/nIrisPerimPoints:2*pi-(2*pi/nIrisPerimPoints);
    irisPoints(1:nIrisPerimPoints,3) = ...
        sin(perimeterPointAngles)*sceneGeometry.eye.irisRadius + sceneGeometry.eye.irisCenter(3);
    irisPoints(1:nIrisPerimPoints,2) = ...
        cos(perimeterPointAngles)*sceneGeometry.eye.irisRadius + sceneGeometry.eye.irisCenter(2);
    irisPoints(1:nIrisPerimPoints,1) = ...
        0 + sceneGeometry.eye.irisCenter(1);
    
    % Add the points and labels
    eyeWorldPoints = [eyeWorldPoints; irisPoints];
    tmpLabels = cell(size(irisPoints,1), 1);
    tmpLabels(:) = {'irisPerimeter'};
    pointLabels = [pointLabels; tmpLabels];
    
    % Create the anterior chamber ellipsoid.
    [p1tmp, p2tmp, p3tmp] = ellipsoid( ...
        sceneGeometry.eye.corneaFrontSurfaceCenter(1), ...
        sceneGeometry.eye.corneaFrontSurfaceCenter(2), ...
        sceneGeometry.eye.corneaFrontSurfaceCenter(3), ...
        sceneGeometry.eye.corneaFrontSurfaceRadius, ...
        sceneGeometry.eye.corneaFrontSurfaceRadius, ...
        sceneGeometry.eye.corneaFrontSurfaceRadius, ...
        30);
    % Convert the surface matrices to a vector of points and switch the
    % axes back
    ansTmp = surf2patch(p1tmp, p2tmp, p3tmp);
    anteriorChamberPoints=ansTmp.vertices;
    
    % Retain those points that are anterior to the iris plane
    retainIdx = logical(...
        (anteriorChamberPoints(:,1) > sceneGeometry.eye.irisCenter(1)));
    if all(~retainIdx)
        error('The pupil plane is set in front of the corneal apea');
    end
    anteriorChamberPoints = anteriorChamberPoints(retainIdx,:);
    
    % Add the points and labels
    eyeWorldPoints = [eyeWorldPoints; anteriorChamberPoints];
    tmpLabels = cell(size(anteriorChamberPoints,1), 1);
    tmpLabels(:) = {'anteriorChamber'};
    pointLabels = [pointLabels; tmpLabels];
    
end

nEyeWorldPoints = size(eyeWorldPoints,1);


%% Project the pupil circle points to headWorld coordinates.
% This coordinate frame is in mm units and has the dimensions (h1,h2,h3).
% The diagram is of a cartoon eye, being viewed directly from the front.
%
%  h1 values negative --> towards the head, positive towards the camera
%
%         h2
%    0,0 ---->
%     |
%  h3 |
%     v
%
%               |
%             __|__
%            /  _  \
%    -------(  (_)  )-------  h2 (horizontal axis of the head)
%            \_____/          rotation about h2 causes pure vertical
%               |             eye movement
%               |
%
%               h3
%   (vertical axis of the head)
%  rotation about h3 causes pure
%     horizontal eye movement
%
%
%
% Position [0,-,-] indicates the front surface of the eye.
% Position [-,0,0] indicates the h2 / h3 position for the center of the
%   pupil when the line that connects the center of the eye and the center
%   of the pupil is normal to the image plane.
%
% We will convert from this coordinate frame to that of the camera scene
% later.


%% Define the eye rotation matrix
% Set pupil torsion to zero, as this will not impact the imaged ellipse
eyeTorsion = 0;

% Assemble a rotation matrix from the head-fixed Euler angle rotations. In
% the head-centered world coordinate frame, positive azimuth, elevation and
% torsion values correspond to leftward, downward and clockwise (as seen
% from the perspective of the subject) eye movements
R3 = [cosd(eyeAzimuth) -sind(eyeAzimuth) 0; sind(eyeAzimuth) cosd(eyeAzimuth) 0; 0 0 1];
R2 = [cosd(eyeElevation) 0 sind(eyeElevation); 0 1 0; -sind(eyeElevation) 0 cosd(eyeElevation)];
R1 = [1 0 0; 0 cosd(eyeTorsion) -sind(eyeTorsion); 0 sind(eyeTorsion) cosd(eyeTorsion)];

% This order (1-2-3) corresponds to a head-fixed, extrinsic, rotation
% matrix. The reverse order (3-2-1) would be an eye-fixed, intrinsic
% rotation matrix and would corresponds to the "Fick coordinate" scheme.
eyeRotation = R1*R2*R3;

%% Prior to rotation, adjust the eyeWorld points for the refraction of the cornea
if isempty(corneaRayTraceFunc)
    virtualEyeWorldPoints = eyeWorldPoints;
else
    % Obtain the functions that calculate the intersection of the ray
    % from an eyeWorldPoint on the camera plane
    [zCameraPlaneX,zCameraPlaneY] = rayIntersectCameraPlane( sceneGeometry, eyeRotation, corneaRayTraceFunc );
    % Find the virtual image points
    virtualEyeWorldPoints = findVirtualImage( eyeWorldPoints, zCameraPlaneX, zCameraPlaneY, corneaRayTraceFunc);
end

% Apply the eye rotation to the pupil plane
headWorldPoints = (eyeRotation*(virtualEyeWorldPoints-sceneGeometry.eye.rotationCenter)')'+sceneGeometry.eye.rotationCenter;


%% Project the pupil circle points to sceneWorld coordinates.
% This coordinate frame is in mm units and has the dimensions (X,Y,Z).
% The diagram is of a cartoon head (borrowed from Leszek Swirski), being
% viewed from above:
%
%   |
%   |    .-.
%   |   |   | <- Head
%   |   `^u^'
% Z |      :V <- Camera    (As seen from above)
%   |      :
%   |      :
%  \|/     o <- Target
%
%     ----------> X
%
% +X = right
% +Y = up
% +Z = front (towards the camera)
%
% The origin [0,0,0] corresponds to the front surface of the eye and the
% center of the pupil when the line that connects the center of rotation of
% the eye and the center of the pupil are normal to the image plane.

% Re-arrange the head world coordinate frame to transform to the scene
% world coordinate frame
sceneWorldPoints = headWorldPoints(:,[2 3 1]);

% We reverse the direction of the Y axis so that positive elevation of the
% eye corresponds to a movement of the pupil upward in the image
sceneWorldPoints(:,2) = sceneWorldPoints(:,2)*(-1);


%% Project the pupil circle points to the image plane
% This coordinate frame is in units of pixels, and has the dimensions
% [x, y]:
%
%      ^
%      |
%   y  |
%      |
%      +------->
% [0,0]    x
%
% With x being left/right and y being down/up
%

% Add a column of ones to support the upcoming matrix multiplication with a
% combined rotation and translation matrix
sceneWorldPoints=[sceneWorldPoints, ones(nEyeWorldPoints,1)];

% Create the projectionMatrix
projectionMatrix = ...
    sceneGeometry.intrinsicCameraMatrix * ...
    [sceneGeometry.extrinsicRotationMatrix, ...
    sceneGeometry.extrinsicTranslationVector];

% Project the world points to the image plane and scale
tmpImagePoints=(projectionMatrix*sceneWorldPoints')';
imagePointsPreDistortion=zeros(nEyeWorldPoints,2);
imagePointsPreDistortion(:,1) = ...
    tmpImagePoints(:,1)./tmpImagePoints(:,3);
imagePointsPreDistortion(:,2) = ...
    tmpImagePoints(:,2)./tmpImagePoints(:,3);


%% Apply radial lens distortion
% This step introduces "pincushion" (or "barrel") distortion produced by
% the lens. The x and y distortion equations are in the normalized image
% coordinates. Thus, the origin is at the optical center (aka principal
% point), and the coordinates are in world units. To apply this distortion
% to our image coordinate points, we subtract the optical center, and then
% divide by fx and fy from the intrinsic matrix.

imagePointsNormalized = (imagePointsPreDistortion - [sceneGeometry.intrinsicCameraMatrix(1,3) sceneGeometry.intrinsicCameraMatrix(2,3)]) ./ ...
    [sceneGeometry.intrinsicCameraMatrix(1,1) sceneGeometry.intrinsicCameraMatrix(2,2)];

% Distortion is proportional to distance from the center of the center of
% projection on the camera sensor
radialPosition = sqrt(imagePointsNormalized(:,1).^2 + imagePointsNormalized(:,2).^2);

distortionVector =   1 + ...
    sceneGeometry.radialDistortionVector(1).*radialPosition.^2 + ...
    sceneGeometry.radialDistortionVector(2).*radialPosition.^4;

imagePointsNormalizedDistorted(:,1) = imagePointsNormalized(:,1).*distortionVector;
imagePointsNormalizedDistorted(:,2) = imagePointsNormalized(:,2).*distortionVector;

imagePoints = (imagePointsNormalizedDistorted .* [sceneGeometry.intrinsicCameraMatrix(1,1) sceneGeometry.intrinsicCameraMatrix(2,2)]) +...
    [sceneGeometry.intrinsicCameraMatrix(1,3) sceneGeometry.intrinsicCameraMatrix(2,3)];


%% Fit the ellipse in the image plane and store values
% Obtain the transparent ellipse params of the projection of the pupil
% circle on the image plane.
pupilEllipseOnImagePlane = ellipse_ex2transparent(...
    ellipse_im2ex(...
    ellipsefit_direct( imagePoints(1:nPupilPerimPoints,1), ...
    imagePoints(1:nPupilPerimPoints,2)  ...
    ) ...
    )...
    );

% place theta within the range of 0 to pi
if pupilEllipseOnImagePlane(5) < 0
    pupilEllipseOnImagePlane(5) = pupilEllipseOnImagePlane(5)+pi;
end

end % pupilProjection_fwd


%% LOCAL FUNCTIONS

function virtualEyeWorldPoints = findVirtualImage( eyeWorldPoints, zCameraPlaneX, zCameraPlaneY, corneaRayTraceFunc )


syms p1 p2 p3
syms theta_p1p2 theta_p1p3

for ii=1:size(eyeWorldPoints,1)
    eyeWorldPoint=eyeWorldPoints(ii,:);
    xErrorFunc = matlabFunction(abs(subs(zCameraPlaneX,[ p1, p2, theta_p1p2],[eyeWorldPoint(1), eyeWorldPoint(2), theta_p1p2])));
    [solution_theta_p1p2, xError] = fminbnd(xErrorFunc,-(deg2rad(45)),(deg2rad(45)));
    yErrorFunc = matlabFunction(abs(subs(zCameraPlaneY,[ p1, p2, p3, theta_p1p2, theta_p1p3],[eyeWorldPoint(1), eyeWorldPoint(2), eyeWorldPoint(3), solution_theta_p1p2, theta_p1p3])));
    [solution_theta_p1p3, yError] = fminbnd(yErrorFunc,-(deg2rad(45)),(deg2rad(45)));
    
    % obtain the
    virtualImageRay = corneaRayTraceFunc(eyeWorldPoint(1), eyeWorldPoint(2), eyeWorldPoint(3), solution_theta_p1p2, solution_theta_p1p3);
    virtualEyeWorldPoints(ii,:) = virtualImageRay(1,:);
end
end


function [zCameraPlaneX,zCameraPlaneY] = rayIntersectCameraPlane( sceneGeometry, eyeRotation, corneaRayTraceFunc )

% Obtain the outputRay within the eye reference frame
syms p1 p2 p3
syms theta_p1p2 theta_p1p3
outputRayEyeWorld = corneaRayTraceFunc(p1,p2,p3,theta_p1p2,theta_p1p3);

% Shift the eyeWorld ray to the rotational center of the eye,
% rotate for this eye pose, undo the centering
outputRayHeadWorld(1,:)=outputRayEyeWorld(1,:)-sceneGeometry.eye.rotationCenter;
outputRayHeadWorld(2,:)=outputRayEyeWorld(2,:)-sceneGeometry.eye.rotationCenter;
outputRayHeadWorld = (eyeRotation*(outputRayHeadWorld)')';
outputRayHeadWorld(1,:)=outputRayHeadWorld(1,:)+sceneGeometry.eye.rotationCenter;
outputRayHeadWorld(2,:)=outputRayHeadWorld(2,:)+sceneGeometry.eye.rotationCenter;

% Re-arrange the head world coordinate frame to transform to the scene
% world coordinate frame
outputRaySceneWorld = outputRayHeadWorld(:,[2 3 1]);

% We reverse the direction of the Y axis so that positive elevation of the
% eye corresponds to a movement of the pupil upward in the image
outputRaySceneWorld(:,2) = outputRaySceneWorld(:,2)*(-1);

% Obtain an expression for X and Y distances between the nodal point of the camera in the sceneWorld plane and the
% point at which the ray will strike the plane that contains the camera
slope_xZ =(outputRaySceneWorld(2,1)-outputRaySceneWorld(1,1))/(outputRaySceneWorld(2,3)-outputRaySceneWorld(1,3));
slope_yZ =(outputRaySceneWorld(2,2)-outputRaySceneWorld(1,2))/(outputRaySceneWorld(2,3)-outputRaySceneWorld(1,3));

zCameraPlaneX = outputRaySceneWorld(1,1)+((sceneGeometry.extrinsicTranslationVector(3)-outputRaySceneWorld(1,3))*slope_xZ);
zCameraPlaneY = outputRaySceneWorld(1,2)+((sceneGeometry.extrinsicTranslationVector(3)-outputRaySceneWorld(1,3))*slope_yZ);

end
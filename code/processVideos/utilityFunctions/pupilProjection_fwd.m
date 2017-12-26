function [imageEllipseParams, pupilCenterImageCoords, pupilCenterSceneCoords] = pupilProjection_fwd(eyeAzimuth, eyeElevation, pupilArea, sceneGeometry)
% Project pupil azimuth, elevation, and area to ellipse on the image plane
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
%   eyeAzimuth            - Horizontal rotation of the eye
%   eyeElevation          - Vertical rotation of the eye
%   pupilArea             - Pupil area (on the eye) in square millimeters.
%                           Not to be confused with the area of the ellipse
%                           on the image plane.
%   sceneGeometry         - A structure that contains the fields:
%                             - eyeRadius: scalar in millimeters
%                             - extrinsicTranslationVector: a 3x1 matrix
%                               in units of millimeters that relates center
%                               of rotation of the eye to the optical axis
%                               of the camera
%                             - extrinsicRotationMatrix: a 3x3 matrix in 
%                               units of millimeters
%                             - intrinsicCameraMatrix: a 3x3 matrix in
%                               arbitrary units (typically pixels) 
%
% Outputs:
%   imageEllipseParams    - Ellipse in transparent form
%   pupilCenterImageCoords  - A 2x1 vector that specifies the x, y location
%                           on the image plane that corresponds to the
%                           center of the pupil
%   pupilCenterSceneCoords - A 3x1 vector that gives the X, Y, Z location
%                           of the center of the pupil in sceneWorld
%                           coordinates. The position [0,0,0] corresponds
%                           to the front surface of the eye and the center
%                           of the pupil when the line that connects the
%                           center of rotation of the eye and the center of
%                           the pupil are normal to the image plane
%

%% Define a pupil circle in pupilWorld coordinates
% This coordinate frame is in mm units and has the dimensions (p1,p2,p3).
% The diagram is of a cartoon pupil, being viewed directly from the front.
%
%  coordinate [0,0,0] corresponds to the point at the center of the pupil.
%  The first dimension is depth, which is unused (set to zero) as we model
%  the points of the pupil as lying on a plane.
%
%  p1 values positive --> closer to the subject / farther away from the
%  camera
%
%                 |
%     ^         __|__
%  +  |        /     \   
% p2  -  -----(   +   )----- 
%  -  |        \_____/      
%     v           |         
%                 |    
%
%            - <-p3-> +
%

% Derive the pupilRadius from the pupilArea
pupilRadius = sqrt(pupilArea./pi);

% We define five points around the pupil circle as five points are needed
% to uniquely define the ellipse in the image plane. A 6th (center) point
% is included so that the discrepancy between the projected center of the
% pupil and measured center of the ellipse may be examined if one is
% curious.
nPerimPoints = 5;
pupilWorldPointsAngles = 0:2*pi/nPerimPoints:2*pi-(2*pi/nPerimPoints);
pupilWorldPoints(:,3) = sin(pupilWorldPointsAngles)*pupilRadius(1);
pupilWorldPoints(:,2) = cos(pupilWorldPointsAngles)*pupilRadius(1);
pupilWorldPoints(:,1) = 0;
pupilWorldPoints(nPerimPoints+1,:) = [0 0 0];


%% Project the pupil circle points to headWorld coordinates.
% This coordinate frame is in mm units and has the dimensions (h1,h2,h3).
% The diagram is of a cartoon eye, being viewed directly from the front.
%
%  h1 values positive --> closer to the subject / farther away from the
%  camera
%
%         h2
%    0,0 ---> 
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

% Define the location of the eye center of rotation in the head-centered
% coordinate frame
centerOfRotation = [-sceneGeometry.eyeRadius 0 0];

% Apply the eye rotation to the pupil plane
headWorldPoints = (eyeRotation*(pupilWorldPoints+centerOfRotation)')'-centerOfRotation;

% We sign reverse the h1 axis (depth) values to conform to our axis
% direction convention
headWorldPoints(:,1)=headWorldPoints(:,1)*(-1);

%% Project the pupil circle points to sceneWorld coordinates.
% This coordinate frame is in mm units and has the dimensions (X,Y,Z).
% The diagram is of a cartoon head (borrowed from Leszek Swirski), being
% viewed from above:
%
%   ^
%   |    .-.
%   |   |   | <- Head
%   |   `^u^'
% Z |      :V <- Camera    (As seen from above)
%   |      :
%   |      :
%   |      o <- Target
%
%     ----------> X
%
% +X = right
% +Y = up
% +Z = back (farther from the camera)
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
sceneWorldPoints=[sceneWorldPoints, ones(nPerimPoints+1,1)];

% Create the projectionMatrix
projectionMatrix = ...
    sceneGeometry.intrinsicCameraMatrix * ...
    [sceneGeometry.extrinsicRotationMatrix, ...
    sceneGeometry.extrinsicTranslationVector];

% Project the world points to the image plane
imagePointsReconstructedUnscaled=(projectionMatrix*sceneWorldPoints')';
imagePointsReconstructed=zeros(nPerimPoints+1,2);
imagePointsReconstructed(:,1) = ...
    imagePointsReconstructedUnscaled(:,1)./imagePointsReconstructedUnscaled(:,3);
imagePointsReconstructed(:,2) = ...
    imagePointsReconstructedUnscaled(:,2)./imagePointsReconstructedUnscaled(:,3);


%% Fit the ellipse in the image plane and store values
% Obtain the transparent ellipse params of the projection of the pupil
% circle on the image plane. 
imageEllipseParams = ellipse_ex2transparent(...
    ellipse_im2ex(...
        ellipsefit_direct(imagePointsReconstructed(1:nPerimPoints,1),imagePointsReconstructed(1:nPerimPoints,2)) ...
        )...
    );

% Store the coordinates of the projection of the center of the pupil on the
% image plane.
pupilCenterImageCoords = ...
    [imagePointsReconstructed(nPerimPoints+1,1) imagePointsReconstructed(nPerimPoints+1,2)];

% Store the coordinates of the center of the pupil in the sceneCoordinate
% space. This is of use for creating Blender movies.
pupilCenterSceneCoords = ...
    [sceneWorldPoints(nPerimPoints+1,1) sceneWorldPoints(nPerimPoints+1,2) sceneWorldPoints(nPerimPoints+1,3)];

end % function


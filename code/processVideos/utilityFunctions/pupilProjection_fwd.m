function [pupilEllipseOnImagePlane, pupilCenterOnImagePlane] = pupilProjection_fwd(eyeParams, sceneGeometry)
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
%                             - eyeRadius: scalar in millimeters
%
% Outputs:
%   pupilEllipseOnImagePlane - A 1x5 vector that contains the parameters of
%                           pupil ellipse on the image plane cast in 
%                           transparent form
%   pupilCenterOnImagePlane - A 1x2 vector that specifies the x, y location
%                           on the image plane corresponding to the center
%                           of the pupil
%


%% Prepare variables
% Separate the eyeParams into individual variables
eyeAzimuth = eyeParams(1);
eyeElevation = eyeParams(2);
pupilRadius = eyeParams(3);


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
%           - <--p3--> +
%

% Five points are defined around the pupil circle, which uniquely
% constrains the ellipse in the image plane. A 6th (center) point is
% included so that the discrepancy between the projected center of the
% pupil and measured center of the ellipse may be examined.
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

% Define the location of the eye center of rotation in the head-centered
% coordinate frame
centerOfRotation = [-sceneGeometry.eyeRadius 0 0];

% Apply the eye rotation to the pupil plane
headWorldPoints = (eyeRotation*(pupilWorldPoints+centerOfRotation)')'-centerOfRotation;

% We sign reverse the h1 axis (depth) values to conform to our axis
% direction convention
headWorldPoints(:,1)=headWorldPoints(:,1)*(-1);

% We sign reverse the h2 and h3 axis values so that azimuth and elevation
% rotations produce the called for directions of movement of the eye.
% (positive azimuth and elevation move the center of the pupil up and to
% the right in the image plane).
headWorldPoints(:,2)=headWorldPoints(:,2)*(-1);
headWorldPoints(:,3)=headWorldPoints(:,3)*(-1);


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
imagePointsUnscaled=(projectionMatrix*sceneWorldPoints')';
imagePoints=zeros(nPerimPoints+1,2);
imagePoints(:,1) = ...
    imagePointsUnscaled(:,1)./imagePointsUnscaled(:,3);
imagePoints(:,2) = ...
    imagePointsUnscaled(:,2)./imagePointsUnscaled(:,3);


%% Apply radial lens distortion
% This step introduces "pincushion" (or "barrel") distortion produced by
% the lens. The x and y distortion equations are in the normalized image
% coordinates. Thus, the origin is at the optical center (aka principal
% point), and the coordinates are in world units. To apply this distortion
% to our image coordinate points, we subtract the optical center, and then
% divide by fx and fy from the intrinsic matrix.

imagePointsNormalized = (imagePoints - [sceneGeometry.intrinsicCameraMatrix(1,3) sceneGeometry.intrinsicCameraMatrix(2,3)]) ./ ...
    [sceneGeometry.intrinsicCameraMatrix(1,1) sceneGeometry.intrinsicCameraMatrix(2,2)];

% Distortion is proportional to distance from the center of the center of
% projection on the camera sensor
radialPosition = sqrt(imagePointsNormalized(:,1).^2 + imagePointsNormalized(:,2).^2); 

distortionVector =   1 + ...
    sceneGeometry.radialDistortionVector(1).*radialPosition.^2 + ...
    sceneGeometry.radialDistortionVector(2).*radialPosition.^4;

imagePointsNormalizedDistorted(:,1) = imagePointsNormalized(:,1).*distortionVector;
imagePointsNormalizedDistorted(:,2) = imagePointsNormalized(:,2).*distortionVector;

imagePointsDistorted = (imagePointsNormalizedDistorted .* [sceneGeometry.intrinsicCameraMatrix(1,1) sceneGeometry.intrinsicCameraMatrix(2,2)]) +...
    [sceneGeometry.intrinsicCameraMatrix(1,3) sceneGeometry.intrinsicCameraMatrix(2,3)];


%% Fit the ellipse in the image plane and store values
% Obtain the transparent ellipse params of the projection of the pupil
% circle on the image plane. 
pupilEllipseOnImagePlane = ellipse_ex2transparent(...
    ellipse_im2ex(...
        ellipsefit_direct( imagePointsDistorted(1:nPerimPoints,1), ...
                           imagePointsDistorted(1:nPerimPoints,2)  ...
                           ) ...
        )...
    );

% place theta within the range of 0 to pi
if pupilEllipseOnImagePlane(5) < 0
    pupilEllipseOnImagePlane(5) = pupilEllipseOnImagePlane(5)+pi;
end

% Store the coordinates of the projection of the center of the pupil on the
% image plane.
pupilCenterOnImagePlane = ...
    [imagePointsDistorted(nPerimPoints+1,1) imagePointsDistorted(nPerimPoints+1,2)];

end % pupilProjection_fwd


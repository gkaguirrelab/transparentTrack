function sceneGeometry = createSceneGeometry(varargin)
% Create and return a sceneGeometry structure
%
% Description:
%   Using default values and passed key/value pairs, this routine creates a
%   sceneGeometry structure, with fields the describe a camera, an eye, and
%   the geometric relationship between them. The fields are:
%
%   intrinsicCameraMatrix - This matrix has the form:
%
%       [fx  s x0]
%       [0  fy y0]
%       [0   0  1]
%
%   where fx and fy are the focal lengths of the camera in the x and y
%   image dimensions, s is the axis skew, and x0, y0 define the principle
%   offset point. For a camera sensor with square pixels, fx = fy. Ideally,
%   skew should be zero. The principle offset point should be in the center
%   of the sensor. Units are traditionally in pixels. These values can be
%   empirically measured for a camera using a calibration approach
%   (https://www.mathworks.com/help/vision/ref/cameramatrix.html). Note
%   that the values in the matrix returned by the Matlab camera calibration
%   routine must be re-arranged to correspond to the X Y image dimension
%   specifications we use here.
%
%   radialDistortionVector - A two element vector of the form:
%
%       [k1 k2]
%
%   that models the radial distortion introduced the lens. This is an
%   empirically measured property of the camera system.
%
%   extrinsicTranslationVector - a vector of the form:
%
%       [x]
%       [y]
%       [z]
%
%   with the values specifying the location (horizontal, vertical, and
%   depth, respectively) of the principle offset point of the camera in mm
%   relative to the scene coordinate system. We define the origin of the
%   scene coordinate system to be x=0, y=0 along the optical axis of the
%   eye, and z=0 to be the apex of the corneal surface.
%
%   extrinsicRotationMatrix - A 3x3 matrix with the fixed values of:
%
%       [1  0  0]
%       [0 -1  0]
%       [0  0 -1]
%
%   These values invert the axes of the scene coordinate system to produce
%   the direction conventions of the image coordinate system. The
%   projection of pupil circles from scene to image is invariant to
%   rotations of the camera matrix, so these valyes should not require
%   adjustment.
%
%   primaryPosition - A 1x3 vector of:
%
%       [eyeAzimuth, eyeElevation, eyeTorsion]
%
%   that specifies the rotation angles (in head fixed axes) for which the
%   eye is in primary position (as defined by Listing's Law). These values
%   may also be used to define the position at which the subject is
%   fixating the origin point of a stimulus array.
%
%   constraintTolerance - A scalar. The inverse projection from ellipse on
%   the image plane to eye params (azimuth, elevation) imposes a constraint
%   on how well the solution must match the shape of the ellipse (defined
%   by ellipse eccentricity and theta) and the area of the ellipse. This
%   constraint is expressed as a proportion of error, relative to either
%   the largest possible area in ellipse shape or an error in area equal to
%   the area of the ellipse itself (i.e., unity). If the constraint is made
%   too stringent, then in an effort to perfectly match the shape of the
%   ellipse, error will increase in matching the position of the ellipse
%   center. It should be noted that even in a noise-free simulation, it is
%   not possible to perfectly match ellipse shape and area while matching
%   ellipse center position, as the shape of the projection of the pupil
%   upon the image plane deviates from perfectly elliptical due to
%   perspective effects. We find that a value in the range 0.01 - 0.03
%   provides an acceptable compromise in empirical data.
%
%   eye -  This is itself a structure that is returned by the function
%   modelEyeParameters(). The parameters define the anatomical properties
%   of the eye, including the size and shape of the anterior and posterior
%   chamber. These parameters are adjusted for the measured spherical
%   refractive error of the subject and (optionally) measured axial length.
%
% Inputs:
%   none
%
% Optional key/value pairs
%  'intrinsicCameraMatrix' - 3x3 matrix
%  'radialDistortionVector' - 1x2 vector of radial distortion parameters
%  'extrinsicRotationMatrix' - 3x3 matrix
%  'extrinsicTranslationVector' - 3x1 vector
%  'primaryPosition'      - 1x3 vector
%  'constraintTolerance'  - Scalar. Range 0-1. Typical value 0.01 - 0.03
%
% Outputs
%	sceneGeometry         - A structure that contains the components of the
%                           projection model.
%


%% input parser
p = inputParser; p.KeepUnmatched = true;

% Optional analysis params
p.addParameter('intrinsicCameraMatrix',[2600.0 0 320; 0 2600 240; 0 0 1],@isnumeric);
p.addParameter('radialDistortionVector',[0 0],@isnumeric);
p.addParameter('extrinsicTranslationVector',[0; 0; 120],@isnumeric);
p.addParameter('extrinsicRotationMatrix',[1 0 0; 0 -1 0; 0 0 -1],@isnumeric);
p.addParameter('primaryPosition',[0 0 0],@isnumeric);
p.addParameter('constraintTolerance',0.02,@isnumeric);

% parse
p.parse(varargin{:})

sceneGeometry.intrinsicCameraMatrix = p.Results.intrinsicCameraMatrix;
sceneGeometry.radialDistortionVector = p.Results.radialDistortionVector;
sceneGeometry.extrinsicTranslationVector = p.Results.extrinsicTranslationVector;
sceneGeometry.extrinsicRotationMatrix = p.Results.extrinsicRotationMatrix;
sceneGeometry.primaryPosition = p.Results.primaryPosition;
sceneGeometry.eye = modelEyeParameters(varargin{:});
sceneGeometry.constraintTolerance = p.Results.constraintTolerance;

end % createSceneGeometry


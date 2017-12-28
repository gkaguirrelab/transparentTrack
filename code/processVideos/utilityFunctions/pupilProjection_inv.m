function [eyeParams, bestMatchEllipseOnImagePlane, constraintViolation] = pupilProjection_inv(pupilEllipseOnImagePlane, sceneGeometry, varargin)
% Project the ellipse on the image plane to the pupil circle in the scene
%
% Description:
%	Given the sceneGeometry and an ellipse on the image plane, this routine
%   attempts to find parameters of the rotation of the eye and pupil sizes
%   that can account for the eccentricity, theta, and area of the ellipse.
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
%   pupilEllipseOnImagePlane - A 1x5 vector that contains the parameters of
%                           pupil ellipse on the image plane cast in
%                           transparent form
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
% Optional key/value pairs:
%  'x0'                   - Starting point of the search for the eyeParams
%  'absoluteEyeAzimuthUB' - The absolute value of the boundary on azimuth
%  'absoluteEyeElevationUB' - The absolute value of the boundary on
%                           elevation
%  'pupilRadiusBounds'    - A 1x2 vector that contains the lower and upper
%                           bound on pupil radius, in mm.
%  'constraintTolerance'  - A scalar value, expressed as a proportion, that
%                           defines the tolerance for violation of the
%                           nonlinear constraint on ellipse eccentricity
%                           and theta and upon ellipse area. The default
%                           value of 0.01 indicates a 1% threshold for
%                           error in matching of (linearized) eccentricty,
%                           theta, or area. Errors larger than this will
%                           result in the constraintViolation flag being
%                           set to 'true'.
%
% Outputs:
%   eyeParams             - A 1x3 vector provides values for [eyeAzimuth,
%                           eyeElevation, pupilRadius]. Azimuth and
%                           elevation are in units of head-centered
%                           (extrinsic) degrees, and pupil radius is in mm.
%   bestMatchEllipseOnImagePlane - A 1x5 vector that contains the
%                           parameters of
%                           pupil ellipse on the image plane cast in
%                           transparent form. This is the output of the
%                           pupilProjection_fwd model for the sceneGeometry
%                           and the eyeParams
%   constraintViolation   - Scalar. Set to 1 if the best solution violated
%                           the ConstraintTolerance for matching ellipse
%                           area, eccentricity, or theta. Otherwise set to
%                           zero.


%% Parse input
p = inputParser;

% Required input
p.addRequired('pupilEllipseOnImagePlane',@isnumeric);
p.addRequired('sceneGeometry',@isstruct);

% Optional params
p.addParameter('x0',[0 0 2],@isnumeric);
p.addParameter('absoluteEyeAzimuthUB',35,@isnumeric);
p.addParameter('absoluteEyeElevationUB',25,@isnumeric);
p.addParameter('pupilRadiusBounds',[0.5,5],@isnumeric);
p.addParameter('constraintTolerance',0.01,@isnumeric);

% Parse and check the parameters
p.parse(pupilEllipseOnImagePlane, sceneGeometry, varargin{:});


%% Assemble bounds
% Because ellipses are symmetric about their axes, a given eccentricity and
% theta of an ellipse is consistent with two possible eye rotation
% solutions. We finesse this ambiguity by bounding the possible azimuth and
% elevation solutions depending upon the quadrant in which the center of
% the ellipse falls, relative to the center of projection derived from the
% scene geometry.

% Identify the center of projection.
projectionMatrix = ...
    sceneGeometry.intrinsicCameraMatrix * ...
    [sceneGeometry.extrinsicRotationMatrix, ...
    sceneGeometry.extrinsicTranslationVector];

CoP = projectionMatrix*[0 0 0 1]';
CoP(1:2)=CoP(1:2)./CoP(3);
CoP=CoP(1:2);

% Set the bounds on the eyeParams based upon the quadrant of the ellipse
% center
if pupilEllipseOnImagePlane(1) < CoP(1)
    eyeParamsLB(1) = -p.Results.absoluteEyeAzimuthUB; eyeParamsUB(1) = 0;
else
    eyeParamsLB(1) = 0; eyeParamsUB(1) = p.Results.absoluteEyeAzimuthUB;
end
if pupilEllipseOnImagePlane(2) > CoP(2)
    eyeParamsLB(2) = -p.Results.absoluteEyeElevationUB; eyeParamsUB(2) = 0;
else
    eyeParamsLB(2) = 0; eyeParamsUB(2) = p.Results.absoluteEyeElevationUB;
end

% Add the pupil radius constraint
eyeParamsLB(3)=p.Results.pupilRadiusBounds(1);
eyeParamsUB(3)=p.Results.pupilRadiusBounds(2);


%% Perform the search
% We use nested functions for the objective and constraint so that the
% forward pupil projection is only computed once for each iteration of the
% fmincon solver

% Define search options
options = optimoptions(@fmincon,...
    'Display','off', ...
    'ConstraintTolerance',p.Results.constraintTolerance);

% Define variables used in the nested functions
xLast = []; % Last place pupilProjection_fwd was called
nestedTargetEllipse = pupilEllipseOnImagePlane; % the target ellipse params
nestedCandidateEllipse = []; % holds pupilProjection_fwd result at xLast
nestedSceneGeometry = sceneGeometry; % a copy of sceneGeometry

% Define anonymous functions for the objective and constraint
objectiveFun = @objfun; % the objective function, nested below
constraintFun = @constr; % the constraint function, nested below

% Call fmincon
[eyeParams, areaError, ~, output] = ...
    fmincon(objectiveFun, p.Results.x0, [], [], [], [], eyeParamsLB, eyeParamsUB, constraintFun, options);

    function fval = objfun(x)
        if ~isequal(x,xLast) % Check if computation is necessary
            nestedCandidateEllipse = pupilProjection_fwd(x, nestedSceneGeometry);
            xLast = x;
        end
        % Now compute objective function
        fval = abs(nestedTargetEllipse(3)-nestedCandidateEllipse(3));
    end

    function [c,ceq] = constr(x)
        if ~isequal(x,xLast) % Check if computation is necessary
            nestedCandidateEllipse = pupilProjection_fwd(x, nestedSceneGeometry);
            xLast = x;
        end
        % We implement here a constraint upon the eccentricity and theta of
        % the ellipse. The theta and eccentricity of an ellipse can be
        % described as a point in polar coordinates. We express the
        % constraint as the vector distance between these points. Direct
        % minimization of differences in theta is a poor constraint, as
        % differences in theta have reduced meaning at small
        % eccentricities. Because ellipses are symmetric, theta spans the
        % range of 0:pi. Therefore, the theta value is doubled prior to
        % conversion to Cartesian coordinates so that the space wraps at
        % the 0 - pi transition point. Eccentricity is linearized to place
        % the parameter in the range of 0:1. The ceq value is divided by 2,
        % so that the largest possible error is unity.
        
        thetaT = nestedTargetEllipse(5)*2;
        thetaC = nestedCandidateEllipse(5)*2;
        rhoT = 1-sqrt(1-nestedTargetEllipse(4)^2);
        rhoC = 1-sqrt(1-nestedCandidateEllipse(4)^2);
                
        ceq = sqrt(rhoT^2 + rhoC^2 - 2*rhoT*rhoC*cos(thetaT-thetaC))/2;

        % c is unused
        c = [];
    end

% Place the final forward model output of the search into a variable to
% return
bestMatchEllipseOnImagePlane = nestedCandidateEllipse;

% We set the constraintViolation flag to true if either the non-linear
% constraint or the areaError (as a proportion of target ellipse area)
% exceeded the tolerance value.
constraintViolation = double(...
    output.constrviolation > p.Results.constraintTolerance || ...
    (areaError/pupilEllipseOnImagePlane(3)) > p.Results.constraintTolerance ...
    );

end % function -- pupilProjection_inv


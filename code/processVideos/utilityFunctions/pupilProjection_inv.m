function [eyeParams, bestMatchEllipseOnImagePlane, centerError, shapeError, areaError] = pupilProjection_inv(pupilEllipseOnImagePlane, sceneGeometry, varargin)
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
%                             - constraintTolearance: A scalar value, 
%                               expressed as a proportion, that defines the
%                               tolerance for violation of the nonlinear
%                               constraints on ellipse shape and area.
%
% Optional key/value pairs:
%  'x0'                   - Starting point of the search for the eyeParams.
%                           If not defined, the starting point will be
%                           estimated from the coordinates of the ellipse
%                           center.
%  'absoluteEyeAzimuthUB' - The absolute value of the boundary on azimuth
%  'absoluteEyeElevationUB' - The absolute value of the boundary on
%                           elevation
%  'pupilRadiusBounds'    - A 1x2 vector that contains the lower and upper
%                           bound on pupil radius, in mm.
%  'constraintTolerance'  - If passed, this value will over-ride the
%                           value in the sceneGeometry structure. 
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
%   centerError           - Scalar. The Euclidean distance (in pixels)
%                           between the [x, y] center of the
%                           pupilEllipseOnImagePlane and the center of the
%                           bestMatchEllipseOnImagePlane.
%   shapeError            - Scalar. The proportion of error in fitting
%                           ellipse shape, range 0-1.
%   areaError             - Scalar. The proportion of error in fitting
%                           ellipse area; unbounded around zero.


%% Parse input
p = inputParser;

% Required input
p.addRequired('pupilEllipseOnImagePlane',@isnumeric);
p.addRequired('sceneGeometry',@isstruct);

% Optional params
p.addParameter('x0',[],@(x)(isempty(x) | isnumeric(x)));
p.addParameter('absoluteEyeAzimuthUB',35,@isnumeric);
p.addParameter('absoluteEyeElevationUB',25,@isnumeric);
p.addParameter('pupilRadiusBounds',[0.5,5],@isnumeric);
p.addParameter('constraintTolerance',[],@(x)(isempty(x) | isnumeric(x)));

% Parse and check the parameters
p.parse(pupilEllipseOnImagePlane, sceneGeometry, varargin{:});

%% Check inputs and handle immediate exits
if isempty(pupilEllipseOnImagePlane)
    centerError=NaN;
    return
end

%% Assemble bounds and x0
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
    eyeParamsLB(1) = -p.Results.absoluteEyeAzimuthUB; eyeParamsUB(1) = .5;
else
    eyeParamsLB(1) = -.5; eyeParamsUB(1) = p.Results.absoluteEyeAzimuthUB;
end
if pupilEllipseOnImagePlane(2) > CoP(2)
    eyeParamsLB(2) = -p.Results.absoluteEyeElevationUB; eyeParamsUB(2) = .5;
else
    eyeParamsLB(2) = -.5; eyeParamsUB(2) = p.Results.absoluteEyeElevationUB;
end

% Add the pupil radius constraint
eyeParamsLB(3)=p.Results.pupilRadiusBounds(1);
eyeParamsUB(3)=p.Results.pupilRadiusBounds(2);

% If x0 is undefined, we will make a guess
if isempty(p.Results.x0)
    % Probe the forward model to determine how many pixels of change in the
    % location of the pupil ellipse correspond to one degree of rotation.
    probeEllipse=pupilProjection_fwd([1 0 2],sceneGeometry);
    pixelsPerDeg = probeEllipse(1)-CoP(1);
    
    % Estimate the eye azimuth and elevation by the X and Y displacement of
    % the ellipse center from the center of projection. Set the initial
    % guess for the pupil radius to 2 mm.
    x0(1) = (pupilEllipseOnImagePlane(1) - CoP(1))/pixelsPerDeg;
    x0(2) = (CoP(2) - pupilEllipseOnImagePlane(2))/pixelsPerDeg;
    
    % Estimate the pupil radius in pixels, accounting for the eccentricity
    % of the ellipse in the image plane
    ellipseAspectRatio = sqrt(1 - (pupilEllipseOnImagePlane(4)^2));
    pupilRadiusPixels = sqrt(pupilEllipseOnImagePlane(3) / (pi * ellipseAspectRatio));
    
    % Probe the forward model at the estimated Azimuth and Elevation to
    % estimate the pupil radius.
    probeEllipse=pupilProjection_fwd([x0(1) x0(2) 2],sceneGeometry);
    pixelsPerMM = sqrt(probeEllipse(3)/pi)/2;
    
    % Set the initial value for pupil radius in mm
    x0(3) = pupilRadiusPixels/pixelsPerMM;
    
    % Ensure that x0 lies within the bounds with a bit of headroom so that
    % the solver does not get stuck up against a bound.
    boundHeadroom = (eyeParamsUB - eyeParamsLB)*0.05;
    x0=min([eyeParamsUB-boundHeadroom; x0]);
    x0=max([eyeParamsLB+boundHeadroom; x0]);
else
    x0 = p.Results.x0;
end

%% Perform the search
% We use nested functions for the objective and constraint so that the
% forward pupil projection is only computed once for each iteration of the
% fmincon solver

% Define search options
options = optimoptions(@fmincon,...
    'Display','off', ...
    'Algorithm','sqp',...
    'OutputFcn',@outfun, ...
    'ConstraintTolerance',p.Results.constraintTolerance);

% Define variables used in the nested functions
xLast = []; % Last place pupilProjection_fwd was called
nestedTargetEllipse = pupilEllipseOnImagePlane; % the target ellipse params
nestedCandidateEllipse = []; % holds pupilProjection_fwd result at xLast
nestedSceneGeometry = sceneGeometry; % a copy of sceneGeometry
history.x = []; % used by outfun
history.fval = []; % used by outfun
history.constraint = []; % used by outfun

% Obtain the constraintTolerance
if isempty(p.Results.constraintTolerance)
    constraintTolerance = sceneGeometry.constraintTolerance;
else
    constraintTolerance = p.Results.constraintTolerance;
end

% Define anonymous functions for the objective and constraint
objectiveFun = @objfun; % the objective function, nested below
constraintFun = @constr; % the constraint function, nested below

% Call fmincon
[eyeParams, centerError, ~, output] = ...
    fmincon(objectiveFun, x0, [], [], [], [], eyeParamsLB, eyeParamsUB, constraintFun, options);

    function fval = objfun(x)
        if ~isequal(x,xLast) % Check if computation is necessary
            nestedCandidateEllipse = pupilProjection_fwd(x, nestedSceneGeometry);
            xLast = x;
        end
        % Now compute objective function in terms of area error squared
        fval = sqrt((nestedTargetEllipse(1) - nestedCandidateEllipse(1))^2 + ...
            (nestedTargetEllipse(2) - nestedCandidateEllipse(2))^2);
    end

    function [c,ceq] = constr(x)
        if ~isequal(x,xLast) % Check if computation is necessary
            nestedCandidateEllipse = pupilProjection_fwd(x, nestedSceneGeometry);
            xLast = x;
        end
        % c:
        % The theta and eccentricity of an ellipse can be
        % described as a point in polar coordinates. We express the
        % constraint as the vector distance between these points. Direct
        % minimization of differences in theta is a poor constraint, as
        % differences in theta have reduced meaning at small
        % eccentricities. Because ellipses are symmetric, theta spans the
        % range of 0:pi. Therefore, the theta value is doubled prior to
        % conversion to Cartesian coordinates so that the space wraps at
        % the 0 - pi transition point. Eccentricity has a value ranging
        % from zero (circular) to 1 (a fully flattened ellipse). The ceq
        % value is divided by 2, so that the largest possible error is
        % unity.
        
        thetaT = nestedTargetEllipse(5)*2;
        thetaC = nestedCandidateEllipse(5)*2;
        rhoT = 1-sqrt(1-nestedTargetEllipse(4)^2);
        rhoC = 1-sqrt(1-nestedCandidateEllipse(4)^2);
        
        c = sqrt(rhoT^2 + rhoC^2 - 2*rhoT*rhoC*cos(thetaT-thetaC))/2;

        % ceq:
        % proportional difference in ellipse areas
        ceq = (nestedTargetEllipse(3) - nestedCandidateEllipse(3))/nestedTargetEllipse(3);
    end

    function stop = outfun(x,optimValues,state)
        stop = false;
        
        switch state
            case 'init'
                % We could set up a figure here if we wished
            case 'iter'
                % Concatenate current point and objective function
                % value with history. x must be a row vector.
                history.fval = [history.fval; optimValues.fval];
                history.constraint = [history.constraint; optimValues.constrviolation];
                history.x = [history.x; x];
                % Test if we are done
                if (optimValues.fval / (nestedTargetEllipse(3)^2)) < constraintTolerance && ...
                    optimValues.constrviolation < constraintTolerance
                    stop = true;
                end
            case 'done'
            otherwise
        end
    end

%% Now optimize for pupil radius


% Store the params of the best fitting ellipse 
bestMatchEllipseOnImagePlane = nestedCandidateEllipse;

% Store the errors
shapeError = output.constrviolation;
areaError = (pupilEllipseOnImagePlane(3) - bestMatchEllipseOnImagePlane(3))/pupilEllipseOnImagePlane(3);

end % function -- pupilProjection_inv


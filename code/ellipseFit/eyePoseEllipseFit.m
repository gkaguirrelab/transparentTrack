function [eyePose, RMSE] = eyePoseEllipseFit(Xp, Yp, sceneGeometry, rayTraceFuncs, varargin)
% Fit an image plane ellipse by perspective projection of a pupil circle
%
% Syntax:
%  [eyePose, RMSE] = eyePoseEllipseFit(Xp, Yp, sceneGeometry, rayTraceFuncs)
%
% Description:
%   The routine fits points on the image plane based upon the eye
%   parameters (azimuth, elevation, pupil radius) that would produce the
%   best fitting ellipse projected according to sceneGeometry.
%
%   The search is constrained by the upper and lower bounds of the eyePose.
%   The default values specified here represent the physical boundaries of
%   the rotation model. Tighter, biologically informed constraints may be
%   passed by the calling function.
%
% Inputs:
%   Xp, Yp                - Vector of points to be fit
%
% Optional key/value pairs:
%  'x0'                   - Initial guess for the eyePose. The initial
%                           azimuth and elevation is slightly different
%                           from zero, as the routines can become stuck in
%                           local minima for rotation values exactly at
%                           zero.
%  'eyePoseLB'            - Lower bound on the eyePose
%  'eyePoseUB'            - Upper bound on the eyePose
%
% Outputs:
%   eyePose               - A 1x4 matrix containing the best fitting eye
%                           parameters (azimuth, elevation, torsion, pupil
%                           radius)
%   RMSE                  - Root mean squared error of the distance of
%                           boundary point in the image to the fitted
%                           ellipse
%


%% Parse input
p = inputParser;

% Required
p.addRequired('Xp',@isnumeric);
p.addRequired('Yp',@isnumeric);
p.addRequired('sceneGeometry',@isstruct);
p.addRequired('rayTraceFuncs',@(x)(isempty(x) | isstruct(x)));

% Optional
p.addParameter('x0',[],@(x)(isempty(x) | isnumeric(x)));
p.addParameter('eyePoseLB',[-89,-89,0,0.1],@isnumeric);
p.addParameter('eyePoseUB',[89,89,0,4],@isnumeric);

% Parse and check the parameters
p.parse(Xp, Yp, sceneGeometry, rayTraceFuncs, varargin{:});


%% Set bounds and x0
% Identify the center of projection.
eyePoseLB = p.Results.eyePoseLB;
eyePoseUB = p.Results.eyePoseUB;

projectionMatrix = ...
    sceneGeometry.intrinsicCameraMatrix * ...
    [sceneGeometry.extrinsicRotationMatrix, ...
    sceneGeometry.extrinsicTranslationVector];

CoP = projectionMatrix*[0 0 0 1]';
CoP(1:2)=CoP(1:2)./CoP(3);
CoP=CoP(1:2);

meanXp = mean(Xp);
meanYp = mean(Yp);

% Set the bounds on the eyePose based upon the quadrant of the ellipse
% center. We provide half a degree of wiggle in the fit around zero.
if meanXp < CoP(1)
    eyePoseUB(1) = .5;
else
    eyePoseLB(1) = -.5;
end
if meanYp > CoP(2)
    eyePoseUB(2) = .5;
else
    eyePoseLB(2) = -.5;
end

% If x0 is undefined, we make a guess based upon the location of the center
% of the points to be fit
if isempty(p.Results.x0)
    % Probe the forward model to determine how many pixels of change in the
    % location of the pupil ellipse correspond to one degree of rotation.
    % Omit ray-tracing to save time as it has minimal effect upon the
    % position of the center of the ellipse.
    probeEllipse=pupilProjection_fwd([1 0 0 2],sceneGeometry, []);
    pixelsPerDeg = probeEllipse(1)-CoP(1);
    
    % Estimate the eye azimuth and elevation by the X and Y displacement of
    % the ellipse center from the center of projection. Torsion is set to
    % zero
    x0(1) = ((meanXp - CoP(1))/pixelsPerDeg);
    x0(2) = ((CoP(2) - meanYp)/pixelsPerDeg);
    x0(3) = 0;
        
    % Force the angles within bounds
    x0=min([eyePoseUB(1:3); x0]);
    x0=max([eyePoseLB(1:3); x0]);    
    
    % Estimate the pupil radius in pixels
    pupilRadiusPixels = max([abs(max(Xp)-min(Xp)) abs(max(Yp)-min(Yp))])/2;
    
    % Probe the forward model at the estimated pose angles to
    % estimate the pupil radius. Here we do need ray tracing as it
    % has a substantial influence upon the area of the ellipse.
    probeEllipse=pupilProjection_fwd([x0(1) x0(2) x0(3) 2], sceneGeometry, rayTraceFuncs);
    pixelsPerMM = sqrt(probeEllipse(3)/pi)/2;
    
    % Set the initial value for pupil radius in mm
    x0(4) = pupilRadiusPixels/pixelsPerMM;
    
    % If the absolute value of an estimated angle is less than 2 degrees,
    % set the value to close to zero. This is done as fmincon seems to
    % avoid solutions exactly at zero, and this kludge fixes that behavior.
    x0(abs(x0(1:3))<2) = 1e-6;
    
    % Ensure that x0 lies within the bounds with a bit of headroom so that
    % the solver does not get stuck up against a bound.
    boundHeadroom = (eyePoseUB - eyePoseLB)*0.001;
    x0=min([eyePoseUB-boundHeadroom; x0]);
    x0=max([eyePoseLB+boundHeadroom; x0]);
else
    x0 = p.Results.x0;
end


% Define an anonymous function for the objective
myObj = @(x) objfun(x, Xp, Yp, sceneGeometry, rayTraceFuncs);

% define some search options
options = optimoptions(@fmincon,...
    'Display','off');

% Perform the non-linear search .We sometimes obtain a singular matrix
% warning here; turn it off temporarily
warningState = warning;
warning('off','MATLAB:nearlySingularMatrix');
warning('off','MATLAB:singularMatrix');

[eyePose, RMSE, exitFlag] = ...
    fmincon(myObj, x0, [], [], [], [], eyePoseLB, eyePoseUB, [], options);

% If exitFlag==2, we might be in a local minimum; try again starting from
% a position close to the point found by the prior search
if exitFlag == 2
    x0 = eyePose+[1e-6 1e-6 0 1e-6];
    [eyePose, RMSE] = ...
        fmincon(myObj, x0, [], [], [], [], p.Results.eyePoseLB, p.Results.eyePoseUB, [], options);
end

% Restore the warning state
warning(warningState);


end % eyeParamEllipseFit


%% LOCAL FUNCTIONS
function fVal = objfun(x, Xp,Yp, sceneGeometry, rayTraceFuncs)
% Define the objective function
explicitEllipse = ellipse_transparent2ex(pupilProjection_fwd(x, sceneGeometry, rayTraceFuncs));
% This is the RMSE of the distance values of the boundary points to
% the ellipse fit. We check for the case in which the
% explicitEllipse contains NAN values, which can happen when the
% eye pose is such that the border of the pupil would not be
% visible through the cornea. In this case, we return a realMax
% value for the fVal.
if any(isnan(explicitEllipse))
    fVal = realmax;
else
    fVal = sqrt(nanmean(ellipsefit_distance(Xp,Yp,explicitEllipse).^2));
end
end % local objective function


% Taken from the non-linear ellipse fitting routine found within the
% "quadfit" matlab central toolbox

function [d,ddp] = ellipsefit_distance(x,y,p)
% Distance of points to ellipse defined with explicit parameters (center,
% axes and tilt).
% P = b^2*((x-cx)*cos(theta)-(y-cy)*sin(theta))
%   + a^2*((x-cx)*sin(theta)+(y-cy)*cos(theta))
%   - a^2*b^2 = 0

pcl = num2cell(p);
[cx,cy,ap,bp,theta] = pcl{:};

% get foot points
if ap > bp
    a = ap;
    b = bp;
else
    a = bp;
    b = ap;
    theta = theta + 0.5*pi;
    if theta > 0.5*pi
        theta = theta - pi;
    end
end
[xf,yf] = ellipsefit_foot(x,y,cx,cy,a,b,theta);

% calculate distance from foot points
d = realsqrt((x-xf).^2 + (y-yf).^2);

% use ellipse equation P = 0 to determine if point inside ellipse (P < 0)
f = b^2.*((x-cx).*cos(theta)-(y-cy).*sin(theta)).^2 + a^2.*((x-cx).*sin(theta)+(y-cy).*cos(theta)).^2 - a^2.*b^2 < 0;

% convert to signed distance, d < 0 inside ellipse
d(f) = -d(f);

if nargout > 1  % FIXME derivatives
    x = xf;
    y = yf;
    % Jacobian matrix, i.e. derivatives w.r.t. parameters
    dPdp = [ ...  % Jacobian J is m-by-n, where m = numel(x) and n = numel(p) = 5
        b.^2.*cos(theta).*(cos(theta).*(cx-x)-sin(theta).*(cy-y)).*2.0+a.^2.*sin(theta).*(cos(theta).*(cy-y)+sin(theta).*(cx-x)).*2.0 ...
        a.^2.*cos(theta).*(cos(theta).*(cy-y)+sin(theta).*(cx-x)).*2.0-b.^2.*sin(theta).*(cos(theta).*(cx-x)-sin(theta).*(cy-y)).*2.0 ...
        a.*(cos(theta).*(cy-y)+sin(theta).*(cx-x)).^2.*2.0-a.*b.^2.*2.0 ...
        b.*(cos(theta).*(cx-x)-sin(theta).*(cy-y)).^2.*2.0-a.^2.*b.*2.0 ...
        a.^2.*(cos(theta).*(cy-y)+sin(theta).*(cx-x)).*(cos(theta).*(cx-x)-sin(theta).*(cy-y)).*2.0-b.^2.*(cos(theta).*(cy-y)+sin(theta).*(cx-x)).*(cos(theta).*(cx-x)-sin(theta).*(cy-y)).*2.0 ...
        ];
    dPdx = b.^2.*cos(theta).*(cos(theta).*(cx-x)-sin(theta).*(cy-y)).*-2.0-a.^2.*sin(theta).*(cos(theta).*(cy-y)+sin(theta).*(cx-x)).*2.0;
    dPdy = a.^2.*cos(theta).*(cos(theta).*(cy-y)+sin(theta).*(cx-x)).*-2.0+b.^2.*sin(theta).*(cos(theta).*(cx-x)-sin(theta).*(cy-y)).*2.0;
    
    % derivative of distance to foot point w.r.t. parameters
    ddp = bsxfun(@rdivide, dPdp, realsqrt(dPdx.^2 + dPdy.^2));
end
end % ellipsefit_distance


function [xf,yf] = ellipsefit_foot(x,y,cx,cy,a,b,theta)
% Foot points obtained by projecting a set of coordinates onto an ellipse.

xfyf = quad2dproj([x y], [cx cy], [a b], theta);
xf = xfyf(:,1);
yf = xfyf(:,2);
end % ellipsefit_foot


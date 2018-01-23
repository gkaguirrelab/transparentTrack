function [eyeParams, RMSE] = eyeParamEllipseFit(Xp, Yp, sceneGeometry, rayTraceFuncs, varargin)
% Fit an image plane ellipse by perspective projection of a pupil circle
%
% Description:
%   The routine fits points on the image plane based upon the eye
%   parameters (azimuth, elevation, pupil radius) that would produce the
%   best fitting ellipse projected according to sceneGeometry.
%
%   The search is constrained by the upper and lower bounds of the
%   eyeParams. The default values specified here represent the physical
%   boundaries of the rotation model. Tighter, biologically informed 
%   constraints may be passed by the calling function.
%
% Inputs:
%   Xp, Yp                - Vector of points to be fit
%
% Optional key/value pairs:
%  'x0'                   - Initial guess for the eyeParams
%  'eyeParamsLB'          - Lower bound on the eyeParams
%  'eyeParamsUB'          - Upper bound on the eyeParams
%
% Outputs:
%   eyeParams             - A 1x4 matrix containing the best fitting eye
%                           parameters (azimuth, elevation, torsion, pupil radius)
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

p.addParameter('x0',[0 0 2],@isnumeric);
p.addParameter('eyeParamsLB',[-89,-89,0,0.1],@isnumeric);
p.addParameter('eyeParamsUB',[89,89,0,4],@isnumeric);

% Parse and check the parameters
p.parse(Xp, Yp, sceneGeometry, rayTraceFuncs, varargin{:});


%% Define the objective function
% This is the RMSE of the distance values of the boundary points to the
% ellipse fit
myFun = @(p) ...
    sqrt(...
        nanmean(...
            ellipsefit_distance(...
                Xp,...
                Yp,...
                ellipse_transparent2ex(...
                    pupilProjection_fwd(p, sceneGeometry, rayTraceFuncs)...
            	)...
        	).^2 ...
    	)...
    );


% define some search options
options = optimoptions(@fmincon,...
    'Display','off');

% Perform the non-linear search
[eyeParams, RMSE] = ...
    fmincon(myFun, p.Results.x0, [], [], [], [], p.Results.eyeParamsLB, p.Results.eyeParamsUB, [], options);

end % eyeParamEllipseFit


%% LOCAL FUNCTIONS
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


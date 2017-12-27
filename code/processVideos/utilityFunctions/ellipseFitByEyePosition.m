function [eyeParams, RMSE] = ellipseFitByEyePosition(Xp, Yp, sceneGeometry, varargin)
% Non-linear fitting of an ellipse to a set of points
%
% Description:
%   The routine fits an ellipse to data by minimizing point-to-curve
%   distance, using an iterative procedure. The search is conducted over
%   the "transparent" ellipse parameterization, which has explicit values
%   for area eccentricity (aspect ratio), and theta. This allows us to set
%   boundaries and non-linear constraints upon these aspects of the fit.
%
%   This routine is a heavily modified version of a non-linear ellipse
%   fitting routine found within the "quadfit" matlab central toolbox. This
%   routine is dependent upon the quadfit toolbox.
%
% Input:
%   Xp, Yp    - Vector of points to be fit
%   lb, ub    - Upper and lower bounds for the fit search, in transparent
%               ellipse form
%   nonlinconst - Function handle to a non-linear constraint function. This
%               function should take as input the set of ellipse parameters
%               in transparent form and return [c, ceq], where the
%               optimizer constrains the solution such that c<=0 and ceq=0.
%               This is an optional input or can be sent as empty.
%
% Output:
%   transparentEllipseParams - Parameters of the best fitting ellipse
%               expressed in transparent form [5x1 vector]
%   RMSE      - Root mean squared error of the distance of each point in
%               the data to the fitted ellipse
%   constraintError - The value of the nonlinear constraint function for
%               the best fitting ellipse
%
%


%% Parse input
p = inputParser;

% Required
p.addRequired('Xp',@isnumeric);
p.addRequired('Yp',@isnumeric);
p.addRequired('sceneGeometry',@isstruct);

p.addParameter('x0',[0 0 2],@isnumeric);
p.addParameter('lb',[-35,-25,0.5],@isnumeric);
p.addParameter('ub',[35,25,4],@isnumeric);
%p.addParameter('nonlinconst',@(x) (isempty(x) || isa(x, 'function_handle')) );

% Parse and check the parameters
p.parse(Xp, Yp, sceneGeometry, varargin{:});


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
                    pupilProjection_fwd(...
                        p(1), p(2), p(3), sceneGeometry...
                    )...
            	)...
        	).^2 ...
    	)...
    );


% define some search options
options = optimoptions(@fmincon,...
    'Display','off');
%     'Algorithm','interior-point',...
%     'Diagnostics','off',...
%     'DiffMinChange', 0.001);

% save the current warning status and silence anticipated warnings
warningState = warning;
% warning('off','MATLAB:nearlySingularMatrix');

% Perform the non-linear search
[eyeParams, RMSE] = ...
    fmincon(myFun, p.Results.x0, [], [], [], [], p.Results.lb, p.Results.ub, [], options);

% Restore the warning state
warning(warningState);

end % function -- constrainedEllipseFit



%% LOCAL FUNCTIONS

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


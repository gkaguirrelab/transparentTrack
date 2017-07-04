function [pFitTransparent, pSD, e] = constrainedEllipseFit(x, y, lb, ub, nonlinconst)
% function [pFitTransparent, pSD, e] = constrainedEllipseFit(x, y, lb, ub, nonlinconst) 
%
% This routine is a modification of a non-linear ellipse fitting routine
% that is found within the "quadfit" matlab central toolbox. This routine
% is dependent upon the quadfit toolbox.
%
% The routine fits an ellipse to data by minimizing point-to-curve
% distance, using an iterative procedure. The fitting approach used here
% has several customizations designed to assist the fitting of the boundary
% of the pupil:
%
%  - The search is conducted over a parameterization of the ellipse that we
%  refer to as "transparent" parameters. The transparent parameter set has
%  explicit value for area and eccentricity (aspect ratio). This allows us
%  to set boudnaries and non-linear constraints upon these aspects of the
%  fit.
%  - A standard deviation of the parameters is estimated from the inverse
%  of the Hessian matrix. This is a flawed estimate for several reasons
%  (not least of which is that this is a constrained search). Nonetheless,
%  this value is found to be useful in subsequent Bayesian smoothing
%  approaches that take place outside of this routine.
%  - As the Hessian is sometimes not invertible, the routine detects for a
%  singular or close to singular matrix and then will use a pseudo-inverse.
%
% Output arguments:
% pFitTransparent:
%    parameters of ellipse expressed in transparent form (row vec)
% pSD:
%    standard deviations of the parameters estimated from the Hessian (row vec)
% e:
%    sqrt of the sum squared distance of the likelihood fit to the data
%
% Input arguments
% x,y:
%    vector of points that define the edge of the pupil to be fit
%
% lb, ub:
%    upper and lower bounds for the fit search (in ellipse transparent
%    form)
%
% nonlinconst:
%    Function handle to a non-linear constraint function. This function
%    should take as input the set of ellipse parameters in transparent form
%    and return [c, ceq], where the optimizer constrains the solution such
%    that c<=0 and ceq=0. This is an optional input or can be sent as
%    empty.

% compute a close-enough initial estimate
pInitImplicit = quad2dfit_taubin(x,y);
switch imconic(pInitImplicit,0)
    case 'ellipse'  % use Taubin's fit, which has produced an ellipse
    otherwise  % use direct least squares ellipse fit to obtain an initial estimate
        pInitImplicit = ellipsefit_direct(x,y);
end

% convert the initial estimate from implicit form to transparent form
pInitTransparent = ellipse_ex2transparent(ellipse_im2ex(pInitImplicit));

% Find the ellipse parameters (in transparent form) using a non-linear
% search.

options = optimset('fmincon');
options = optimset(options,'Diagnostics','off','Display','off','LargeScale','off','Algorithm','interior-point');

% Define the objective function, which is the sqrt of the sum of squared
% distance values of the boundary points to the ellipse fit
myFun = @(p) sqrt(nansum(ellipsefit_distance(x,y,ellipse_transparent2ex(p)).^2));

% save the current warning status and silence anticipated warnings
warningState = warning;
warning('off','MATLAB:nearlySingularMatrix');

% Fit that sucker
[pFitTransparent,e,~,~,~,~,Hessian] = fmincon(myFun, pInitTransparent, [], [], [], [], lb, ub, nonlinconst, options);

% Restore the warning state
warning(warningState);

% The sqrt of the diagonals of the inverse Hessian matrix approxiates the
%  standard deviation of the parameter estimates (with multiple caveats
%  regarding how the Hessian returned by fmincon is inaccurate for this
%  purpose).

% We check to see if there was a warning regardng an inability to invert
% the matrix, and if so we used the pseudo-inverse

% save the current warning status and silence anticipated warnings
warningState = warning;
warning('off','MATLAB:nearlySingularMatrix');
warning('off','MATLAB:singularMatrix');

% empty the warning string
lastwarn('');

% try the inverse
pSD = sqrt(diag(inv(Hessian)));

% get the lastwarn status and used pinv if there was a problem
[~, msgid] = lastwarn;
switch msgid
    case 'MATLAB:nearlySingularMatrix'
      pSD = sqrt(diag(pinv(Hessian)));
    case 'MATLAB:singularMatrix'
      pSD = sqrt(diag(pinv(Hessian)));
end

% restore the warning status
warning(warningState);

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

function [d,ddp] = ellipsefit_distance_kepler(x,y,p)
% Distance of points to ellipse defined with Kepler's parameters.

pcl = num2cell(p);
[px,py,qx,qy,a] = pcl{:};

% get foot points
[xf,yf] = ellipsefit_foot_kepler(x,y,px,py,qx,qy,a);

% calculate distance from foot points
d = realsqrt((x-xf).^2 + (y-yf).^2);

% use ellipse equation P = 0 to determine if point inside ellipse (P < 0)
f = realsqrt((x-px).^2+(y-py).^2) + realsqrt((x-qx).^2+(y-qy).^2) - 2*a < 0;

% convert to signed distance, d < 0 inside ellipse
d(f) = -d(f);

if nargout > 1
    % compute partial derivatives of ellipse equation P given with Kepler's parameters
    % in Kepler's form, the foci of an ellipse are (px;py) and (qx;qy),
    % and 2*a is the major axis such that the parameter vector is [p1 p2 q1 q2 a]
    
    % Jacobian matrix, i.e. derivatives w.r.t. parameters
    dPdp = [ ...  % Jacobian J is m-by-n, where m = numel(x) and n = numel(p) = 5
        -(xf-px)./realsqrt((xf-px).^2+(yf-py).^2), ...
        -(yf-py)./realsqrt((xf-px).^2+(yf-py).^2), ...
        -(xf-qx)./realsqrt((xf-qx).^2+(yf-qy).^2), ...
        -(yf-qy)./realsqrt((xf-qx).^2+(yf-qy).^2), ...
        -2.*ones(size(xf)) ...
        ];
    dPdx = (xf-px)./realsqrt((xf-px).^2+(yf-py).^2) + (xf-qx)./realsqrt((xf-qx).^2+(yf-qy).^2);
    dPdy = (yf-py)./realsqrt((xf-px).^2+(yf-py).^2) + (yf-qy)./realsqrt((xf-qx).^2+(yf-qy).^2);
    
    % derivative of distance to foot point w.r.t. parameters
    ddp = bsxfun(@rdivide, dPdp, realsqrt(dPdx.^2 + dPdy.^2));
end

function [xf,yf] = ellipsefit_foot(x,y,cx,cy,a,b,theta)
% Foot points obtained by projecting a set of coordinates onto an ellipse.

xfyf = quad2dproj([x y], [cx cy], [a b], theta);
xf = xfyf(:,1);
yf = xfyf(:,2);

function [xf,yf] = ellipsefit_foot_kepler(x,y,px,py,qx,qy,a)
% Foot points obtained by projecting a set of coordinates onto an ellipse.
%
% Input arguments:
% x,y:
%    coordinates to data points to project
% px,py,qx,qy:
%    coordinates of ellipse foci
% a:
%    ellipse semi-major axis length

c1 = 0.5*(px+qx);
c2 = 0.5*(py+qy);
cf2 = (c1-px)^2 + (c2-py)^2;  % distance squared from center to focus
b = realsqrt(a^2 - cf2);
theta = atan2(qy-py,qx-px);  % tilt angle
[xf,yf] = ellipsefit_foot(x,y,c1,c2,a,b,theta);
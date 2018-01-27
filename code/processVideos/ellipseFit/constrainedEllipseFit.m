function [transparentEllipseParams, RMSE, constraintError] = constrainedEllipseFit(Xp, Yp, lb, ub, nonlinconst)
% Non-linear fitting of an ellipse to a set of points
%
% Syntax:
%  [transparentEllipseParams, RMSE, constraintError] = constrainedEllipseFit(Xp, Yp, lb, ub, nonlinconst)
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
%	Xp, Yp                - Vector of points to be fit
%   lb, ub                - Upper and lower bounds for the fit search, in 
%                           transparent ellipse form
%   nonlinconst           - Function handle to a non-linear constraint
%                           function. This function should take as input
%                           the set of ellipse parameters in transparent
%                           form and return [c, ceq], where the optimizer
%                           constrains the solution such that c<=0 and
%                           ceq=0. This is an optional input or can be sent
%                           as empty.
%
% Output:
%   transparentEllipseParams - Parameters of the best fitting ellipse
%                           expressed in transparent form [5x1 vector]
%   RMSE                  - Root mean squared error of the distance of each
%                           point in the data to the fitted ellipse
%   constraintError       - The value of the nonlinear constraint function
%                           for the best fitting ellipse
%
%


%% Parse input
p = inputParser;

% Required
p.addRequired('Xp',@isnumeric);
p.addRequired('Yp',@isnumeric);
p.addRequired('ub',@isnumeric);
p.addRequired('lb',@isnumeric);
p.addRequired('nonlinconst',@(x) (isempty(x) || isa(x, 'function_handle')) );

% Parse and check the parameters
p.parse(Xp, Yp, ub, lb, nonlinconst);


%% Make an initial guess at the ellipse parameters
% This attempt is placed in a try-catch block, as the attempt can fail and
% return non-real numbers.
try
    % We sometimes obtain a singular matrix warning here; turn it
    % off temporarily
    warningState = warning;
    warning('off','MATLAB:singularMatrix');
    % use direct least squares ellipse fit to obtain an initial
    % estimate
    pInitImplicit = ellipsefit_direct(Xp,Yp);
    % Restore the warning state
    warning(warningState);
    % convert the initial estimate from implicit form to transparent form
    pInitTransparent = ellipse_ex2transparent(ellipse_im2ex(pInitImplicit));
    % place theta within the range of 0 to pi
    if pInitTransparent(5) < 0
        pInitTransparent(5) = pInitTransparent(5)+pi;
    end
catch
    % We couldn't find anything vaguely elliptical; return nans
    transparentEllipseParams=nan(1,5);
    RMSE=nan;
    constraintError=nan;
    return
end


%% Define the objective function
% This is the RMSE of the distance values of the boundary points to the
% ellipse fit
myFun = @(p) sqrt(nanmean(ellipsefit_distance(Xp,Yp,ellipse_transparent2ex(p)).^2));

% If the bounds and the nonlinear constraint function are all empty, then
% just return the initial estimate (and RMSE) obtained by direct fitting
if isempty(ub) && isempty(lb) && isempty(nonlinconst)
    transparentEllipseParams = pInitTransparent;
    RMSE = myFun(transparentEllipseParams);
    constraintError = nan;
    return
end


%% Perform non-linear search for transparent ellipse params

% Force the X and Y values of the initial guess to satisfy the nonlinear
% constraint
if ~isempty(nonlinconst)
    [~, ~, projectedEllipseOnImagePlane] = nonlinconst(pInitTransparent);
    pInitTransparent(4:5)=projectedEllipseOnImagePlane(4:5);
end

% define some search options
options = optimoptions(@fmincon,...
    'Algorithm','interior-point',...
    'Diagnostics','off',...
    'Display','off',...
    'DiffMinChange', 0.001);

% save the current warning status and silence anticipated warnings
warningState = warning;
warning('off','MATLAB:nearlySingularMatrix');

% Perform the non-linear search
[transparentEllipseParams, RMSE, ~, output] = ...
    fmincon(myFun, pInitTransparent, [], [], [], [], lb, ub, nonlinconst, options);

% Extract the constraint error from the output structure
constraintError = output.constrviolation;

% If we are close to zero for eccentricity, we may be in a local minimum
% that tends to find circles. Try a multiStart search with an increased
% lower bound on eccentricity.
if transparentEllipseParams(4) < 1e-12
    adjustedLB = lb;
    adjustedLB(4) = min([max([lb(4) 0.1]) ub(4)]);
    problem = createOptimProblem('fmincon','x0',pInitTransparent,...
        'objective',myFun,'lb',adjustedLB,'ub',ub,...
        'nonlcon',nonlinconst,'options',options);
    ms = MultiStart('Display','off','MaxTime',30,'StartPointsToRun','bounds-ineqs');
    [mstransparentEllipseParams, msRMSE, ~, msOutput] = run(ms,problem,20);
    if msRMSE <= RMSE
        RMSE = msRMSE;
        transparentEllipseParams = mstransparentEllipseParams;
        constraintError = msOutput.constrviolation;
    end
end

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


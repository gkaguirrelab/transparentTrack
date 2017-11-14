function [transparentEllipseParams, RMSE] = constrainedEllipseFit(x, y, lb, ub, nonlinconst, varargin)
% constrainedEllipseFit(x, y, lb, ub, nonlinconst)
%
% This routine is a modification of a non-linear ellipse fitting routine
% that is found within the "quadfit" matlab central toolbox. This routine
% is dependent upon the quadfit toolbox.
%
% The routine fits an ellipse to data by minimizing point-to-curve
% distance, using an iterative procedure. The search is conducted over a
% parameterization of the ellipse that we refer to as "transparent"
% parameters. The transparent parameter set has explicit value for area and
% eccentricity (aspect ratio). This allows us to set boudnaries and
% non-linear constraints upon these aspects of the fit.
%
% Output arguments:
%   transparentEllipseParams - parameters of ellipse expressed in transparent
%       form (row vec)
%   RMSE - root mean squared error of the distance of each
%       point in the data to the fitted ellipse
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


%% Parse input
p = inputParser;

% Required
p.addRequired('x',@isnumeric);
p.addRequired('y',@isnumeric);
p.addRequired('ub',@isnumeric);
p.addRequired('lb',@isnumeric);
p.addRequired('nonlinconst',@(x) (isempty(x) || isa(x, 'function_handle')) );

% Optional analysis params
p.addParameter('globalSearchMaxTimeSeconds', 5, @isnumeric);


% Parse and check the parameters
p.parse(x, y, ub, lb, nonlinconst, varargin{:});


%% Calculate an initial estimate of the ellipse parameters
% This attempt is placed in a try-catch block, as the attempt can fail and
% return non-real numbers.
try
    pInitImplicit = quad2dfit_taubin(x,y);
    switch imconic(pInitImplicit,0)
        case 'ellipse'
            % the initial fit returned an ellipse using Taubin's method. No
            % further action is required.
        otherwise
            % The initial fit with Taubin's method didn't work, so try a
            % direct ellipse fit
            
            % We someties obtain a singular matrix warning here; turn it
            % off temporarily
            warningState = warning;
            warning('off','MATLAB:singularMatrix');
            % use direct least squares ellipse fit to obtain an initial
            % estimate
            pInitImplicit = ellipsefit_direct(x,y);
            % Restore the warning state
            warning(warningState);
    end
    % convert the initial estimate from implicit form to transparent form
    pInitTransparent = ellipse_ex2transparent(ellipse_im2ex(pInitImplicit));
catch
    % We couldn't find anything vaguely elliptical; return nans
    transparentEllipseParams=nan(1,5);
    RMSE=nan;
    return
end


%% Perform non-linear search for transparent ellipse params

% define some search options
options = optimset('fmincon');
options = optimset(options,'Diagnostics','off','Display','off','LargeScale','off','Algorithm','sqp');

% Define the objective function, which is the RMSE of the distance values
% of the boundary points to the ellipse fit
myFun = @(p) sqrt(nanmean(ellipsefit_distance(x,y,ellipse_transparent2ex(p)).^2));

% save the current warning status and silence anticipated warnings
warningState = warning;
warning('off','MATLAB:nearlySingularMatrix');

% Perform the non-linear search
[transparentEllipseParams, RMSE, ~, output] = ...
    fmincon(myFun, pInitTransparent, [], [], [], [], lb, ub, nonlinconst, options);

% If we received a local minimum warning on exit, re-fit the ellipse using
% a more time-consuming global minimization approach.
if ~isempty(strfind(output.message,'Local minimum'))
    problem = createOptimProblem('fmincon','x0',pInitTransparent,...
        'objective',myFun,'lb',lb,'ub',ub,...
        'nonlcon',nonlinconst,'options',options);
    gs = GlobalSearch('Display','off','MaxTime',p.Results.globalSearchMaxTimeSeconds,'StartPointsToRun','bounds-ineqs');
    [transparentEllipseParams, RMSE] = run(gs,problem);
end

% Restore the warning state
warning(warningState);

end % MAIN -- constrainedEllipseFit



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


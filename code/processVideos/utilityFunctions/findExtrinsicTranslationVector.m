function [sceneGeometry, rmseDistanceError] = findExtrinsicTranslationVector(targetEllipses, initialSceneGeometry)

sceneGeometry = initialSceneGeometry;

for pp=1:3
    [extrinsicTranslationVector, eyeRadius, rmseDistanceError]  = performSearch(targetEllipses, sceneGeometry);
    sceneGeometry.extrinsicTranslationVector = extrinsicTranslationVector;
    sceneGeometry.eyeRadius = eyeRadius;
end % for loop

end % function


function [extrinsicTranslationVector, eyeRadius, rmseDistanceError]  = performSearch(targetEllipses, sceneGeometry)

% Setup the initial guess
x0 = [sceneGeometry.extrinsicTranslationVector; sceneGeometry.eyeRadius];
lb = [-5, -5, 20, 11];
ub = [5, 5, 50, 14];
nStarts = 16;

% Define fmincon search options
options = optimoptions(@fmincon,...
    'Display','off', ...
    'UseParallel',true, ...
    'ConstraintTolerance',0.01);

% Define anonymous functions for the objective and constraint
objectiveFun = @objfun; % the objective function, nested below
constraintFun = @constr; % the constraint function, nested below

% Define nested variables for within the search
nEllipses = size(targetEllipses,1);
xLast = []; % Last place pupilProjection_fwd was called
centerErrors = zeros(nEllipses,1);
shapeErrors = zeros(nEllipses,1);

% Set up the search object and problem
searchObj = MultiStart('Display','iter', ...
    'UseParallel',true, ...
    'MaxTime',300);

problem = createOptimProblem('fmincon',...
    'x0',x0, ...
    'objective',objectiveFun, ...
    'lb',lb,'ub',ub,...
    'options',options, ...
    'nonlcon',constraintFun);

[x, rmseDistanceError] = ...
    run(searchObj, problem, nStarts);
    function fval = objfun(x)
        if ~isequal(x,xLast) % Check if computation is necessary
            candidateSceneGeometry = sceneGeometry;
            candidateSceneGeometry.extrinsicTranslationVector = x(1:3);
            candidateSceneGeometry.eyeRadius = x(4);
            for ii=1:nEllipses
                [~, ~, centerErrors(ii), shapeErrors(ii), ~] = ...
                    pupilProjection_inv(targetEllipses(ii,:), candidateSceneGeometry, 'constraintTolerance', 0.01);
            end
            xLast = x;
        end
        
        % Now compute objective function as the RMSE of the distance
        % between the taget and modeled ellipses
        fval = sqrt(mean(centerErrors.^2));
    end

    function [c,ceq] = constr(x)
        if ~isequal(x,xLast) % Check if computation is necessary
            candidateSceneGeometry = sceneGeometry;
            candidateSceneGeometry.extrinsicTranslationVector = x(1:3);
            candidateSceneGeometry.eyeRadius = x(4);
            for ii=1:nEllipses
                [~, ~, centerErrors(ii), shapeErrors(ii), ~] = ...
                    pupilProjection_inv(targetEllipses(ii,:), candidateSceneGeometry, 'constraintTolerance', 0.01);
            end
            xLast = x;
        end
        
        % The constraint is the max of constraint violations across target
        % ellipses
        %        ceq = max(shapeErrors);
        ceq = [];
        
        % c is unused
        c = [];
    end

extrinsicTranslationVector = x(1:3);
eyeRadius = x(4);

end %
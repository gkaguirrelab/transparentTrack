function [bestFitSceneGeometry, distanceError, constraintViolation] = findExtrinsicTranslationVector(targetEllipses, initialSceneGeometry)

% Define search options
options = optimoptions(@fmincon,...
    'Display','iter', ...
    'UseParallel',true, ...
    'ConstraintTolerance',0.01);

nEllipses = size(targetEllipses,1);

xLast = []; % Last place pupilProjection_fwd was called
distanceErrors = zeros(nEllipses,1);
constraintViolations = zeros(nEllipses,1);

% Define anonymous functions for the objective and constraint
objectiveFun = @objfun; % the objective function, nested below
constraintFun = @constr; % the constraint function, nested below

% Call fmincon
[extrinsicTranslationVector, distanceError, ~, output] = ...
    fmincon(objectiveFun, initialSceneGeometry.extrinsicTranslationVector, [], [], [], [], [-5, -5, 20], [5, 5, 50], constraintFun, options);

    function fval = objfun(x)
        if ~isequal(x,xLast) % Check if computation is necessary
            candidateSceneGeometry = initialSceneGeometry;
            candidateSceneGeometry.extrinsicTranslationVector = x(1:3);
            for ii=1:nEllipses
                [~, bestMatchEllipseOnImagePlane, shapeError, areaError] = pupilProjection_inv(targetEllipses(ii,:), candidateSceneGeometry);
                constraintViolations(ii) = shapeError + areaError;
                distanceErrors(ii) = sqrt((targetEllipses(ii,1) - bestMatchEllipseOnImagePlane(1))^2 + ...
                    (targetEllipses(ii,2) - bestMatchEllipseOnImagePlane(2))^2);
            end
            xLast = x;
        end
        
        % Now compute objective function as the RMSE of the distance
        % between the taget and modeled ellipses
        fval = sqrt(mean(distanceErrors.^2));
    end

    function [c,ceq] = constr(x)
        if ~isequal(x,xLast) % Check if computation is necessary
            candidateSceneGeometry = initialSceneGeometry;
            candidateSceneGeometry.extrinsicTranslationVector = x(1:3);
            for ii=1:nEllipses
                [~, bestMatchEllipseOnImagePlane, shapeError, areaError] = pupilProjection_inv(targetEllipses(ii,:), candidateSceneGeometry);
                constraintViolations(ii) = shapeError + areaError;
                distanceErrors(ii) = sqrt((targetEllipses(ii,1) - bestMatchEllipseOnImagePlane(1))^2 + ...
                    (targetEllipses(ii,2) - bestMatchEllipseOnImagePlane(2))^2);
            end
            xLast = x;
        end
        
        % The constraint is the sum of constraint violations across target
        % ellipses
        ceq = max(constraintViolations);
        
        % c is unused
        c = [];
    end

bestFitSceneGeometry = initialSceneGeometry;
bestFitSceneGeometry.extrinsicTranslationVector = extrinsicTranslationVector;

constraintViolation = output.constrviolation;

end % function
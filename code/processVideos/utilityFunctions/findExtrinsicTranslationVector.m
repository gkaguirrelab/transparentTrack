function [sceneGeometry, rmseDistanceError] = findExtrinsicTranslationVector(targetEllipses, initialSceneGeometry)

sceneGeometry = initialSceneGeometry;

% Setup the initial guess
x0 = [sceneGeometry.extrinsicTranslationVector; sceneGeometry.eyeRadius];
lb = [-20, -20, 25, 11];
ub = [20, 20, 50, 15];

% Define search options
options = optimoptions(@patternsearch, ...
    'Display','iter',...
    'AccelerateMesh',false,...
    'FunctionTolerance',0.001);

% Define anonymous functions for the objective and constraint
objectiveFun = @objfun; % the objective function, nested below

% Define nested variables for within the search
nEllipses = size(targetEllipses,1);
xLast = []; % Last place pupilProjection_fwd was called
centerErrors = zeros(nEllipses,1);

%[x,rmseDistanceError] = simulannealbnd(objectiveFun,x0,lb,ub,options);
[x,rmseDistanceError] = patternsearch(objectiveFun,x0,[],[],[],[],lb,ub,[],options);
    function fval = objfun(x)
        if ~isequal(x,xLast) % Check if computation is necessary
            candidateSceneGeometry = sceneGeometry;
            candidateSceneGeometry.extrinsicTranslationVector = x(1:3);
            candidateSceneGeometry.eyeRadius = x(4);
            for ii=1:nEllipses
                [~, ~, centerErrors(ii), ~, ~] = ...
                    pupilProjection_inv(targetEllipses(ii,:), candidateSceneGeometry, 'constraintTolerance', 0.01);
            end
            xLast = x;
        end
        
        % Now compute objective function as the RMSE of the distance
        % between the taget and modeled ellipses
        fval = sqrt(mean(centerErrors.^2));
    end

sceneGeometry.extrinsicTranslationVector = x(1:3);
sceneGeometry.eyeRadius = x(4);

end %
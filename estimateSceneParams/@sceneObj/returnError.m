function fVal = returnError(obj, x, errorArgs)
% Updates the scene and head motion models and then returns the objective
%
% Syntax:
%  fVal = obj.returnError(x)
%
% Description:
%   This is the primary objective function of the search across scene
%   parameters for this scene object.
%
% Inputs:
%   x
%   errorArgs             - Cell array of key-values to be sent to
%                           updateError to optionally over-ride the default
%                           parameters.
%
% Outputs:
%   fVal                  - Scalar. The objective function value.
%

% If the stored x is is not empty, and passed x is the same as the current
% state of the model, just return the current fVal to save on computation
% time.
if ~isempty(obj.x)
    if all(obj.x == x)
        fVal = obj.fVal;
        return
    end
    obj.xLast = obj.x;
else
    obj.xLast = nan(size(x));
end

% Update the object with the new parameters
obj.x = x;

% Update the effect of head motion upon relative camera position
obj.updateHead;

% Update the sceneGeometry
obj.updateScene;

% Store the passed errorArgs
obj.errorArgs = errorArgs;

% Update the error
obj.updateError(obj.errorArgs{:});

% Update fValBest and xBest
if obj.fVal < obj.fValBest
    obj.fValBest = obj.fVal;
    obj.xBest = x;
end

% Update xHist and fValHist
obj.xHist(end+1,:) = x;
obj.fValHist(end+1) = obj.fVal;

% Return the fVal
fVal = obj.fVal;

end





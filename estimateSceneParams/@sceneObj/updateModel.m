function fVal = updateModel(obj, x)

% If the stored x is is not empty, and passed x is the same as the current
% state of the model, just return the current fVal
if ~isempty(obj.x)
    if all(obj.x == x)
        fVal = obj.fVal;
        return
    end
end

% Update the object with the new parameters
obj.x = x;

% Update the relative camera position
obj.updateRelCamPos( x(1:4) );

% Update the sceneGeometry
obj.updateScene( x(5:17) );

% Update the error
obj.updateError(obj.keyVals{:});

% Update fValBest and xBest
if obj.fVal < obj.fValBest
    obj.fValBest = obj.fVal;
    obj.xBest = x;
end

fVal = obj.fVal;

end





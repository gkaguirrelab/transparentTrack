function fVal = multiSceneObjective(x,mySceneObjects,model,strategy,verbose)
% An objective function across multiple scenes for estimateSceneParams
%
% Syntax:
%  fVal = multiSceneObjective(x,mySceneObjects,model,strategy,verbose)
%
% Inputs:
%   x                     - 1xp vector of parameters.
%   sceneObjects          - A cell array of handles to objects of the type
%                           sceneObj.
%   model                 - Struct. Describes parameters of the search.
%   strategy              - Char vector. Directs the function to the field
%                              model.strategy.(strategy)
%   verbose               - Logical.
%
% Outputs:
%   fVal                  - Scalar. The objective function value.
%


% Loop over the scenes. Pass the scene params appropriate to that
% scene. Store the objective values.
nScenes = length(mySceneObjects);
fValScene = nan(1,nScenes);
for ss = 1:nScenes
    
    % These are the parameters from x that are appropriate for this scene
    subX = model.func.subX(x,ss);
    
    % Call the object to calculate the updated error function for these
    % parameters
    fValScene(ss) = mySceneObjects{ss}.returnError(subX);
end

% Take the norm of the scene objectives, using the metric defined in the
% model for this strategy
fVal = norm(fValScene,model.strategy.(strategy).multiSceneNorm);

% BADS can't handle Inf in the objective, so replace with real max
fVal = min([fVal realmax]);

% Apply the regularization to penalize changes in the camera depth
penalty = model.func.penalty(x,model.x0,model.strategy.(strategy).penaltyWeight);
fVal = fVal * penalty;

% Each sceneObject has a multiSceneMeta property that is used to stash
% information about the search progress. First determine if we need to
% update the meta information, which is the case if the meta field is
% currently empty, or if the fVal is better than the stored value.
updateMetaFlag = false;
if isempty(mySceneObjects{1}.multiSceneMeta)
    updateMetaFlag = true;
else
    if fVal < mySceneObjects{1}.multiSceneMeta.fVal
        updateMetaFlag = true;
    end
end

% Update the meta value fields if appropriate
if updateMetaFlag    
    
    % Set up the meta data
    multiSceneMeta.x = x;
    multiSceneMeta.fValScene = fValScene;
    multiSceneMeta.fVal = fVal;
    multiSceneMeta.penalty = penalty;
    
    % Loop over the scene objects and store the meta data
    for ss = 1:nScenes
        mySceneObjects{ss}.multiSceneMeta = multiSceneMeta;
        mySceneObjects{ss}.multiSceneIdx = ss;
    end
    
    % If we have the verbose flag, report the new, best fVal and params
    if verbose
        str = 'fVal = %2.2d; x = [ ';
        for pp = 1:length(x)-1
            str = [str '%2.2f, '];
        end
        str = sprintf([str '%2.2f ]\n'],fVal,x);
        fprintf(str);
    end
    
end


end
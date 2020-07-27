function fVal = multiSceneObjective(x,sceneObjects,model,strategy,stage,errorArgsIn,verbose)
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
%   stage                 - Scalar. The stage of the search
%   errorArgsIn           - Cell array. Passed errorArgs from the calling 
%                           function.
%   verbose               - Logical.
%
% Outputs:
%   fVal                  - Scalar. The objective function value.
%


% Get the errorArgs for this stage of the search
errorArgs = [model.strategy.(strategy).errorArgs{stage}, errorArgsIn];

% Loop over the scenes.
nScenes = length(sceneObjects);
fValScene = nan(1,nScenes);
for ss = 1:nScenes
    
    % These are the parameters from x that are appropriate for this scene
    subX = model.func.subX(x,ss);
        
    % Call the object to calculate the updated error function for these
    % parameters
    fValScene(ss) = sceneObjects{ss}.returnError(subX, errorArgs);
    
end

% Take the norm of the scene objectives, using the metric defined in the
% model for this strategy
fVal = norm(fValScene,model.strategy.(strategy).multiSceneNorm);

% BADS can't handle Inf in the objective, so replace with real max
fVal = min([fVal realmax]);

% Apply the regularization to penalize changes in the camera depth and
% torsion
penalty = model.func.penalty(x,model.x0,model.strategy.(strategy).penaltyWeight);
fVal = fVal * penalty;

% Each sceneObject has a multiSceneMeta property that is used to stash
% information about the search progress. First determine if we need to
% update the meta information, which is the case if the meta field is
% currently empty, if the fVal is better than the stored value, or if we
% have started a new stage
updateMetaFlag = false;
if isempty(sceneObjects{1}.multiSceneMeta)
    updateMetaFlag = true;
else
    if fVal < sceneObjects{1}.multiSceneMeta.fVal
        updateMetaFlag = true;
    end
    if stage ~= sceneObjects{1}.multiSceneMeta.currentStage
        updateMetaFlag = true;
    end
end

% Update the meta value fields if appropriate
if updateMetaFlag    
    
    % Loop over the scene objects and store the meta data
    for ss = 1:nScenes
        
        % Add the meta data
        sceneObjects{ss}.multiSceneMeta.stages{stage}.x = x;
        sceneObjects{ss}.multiSceneMeta.stages{stage}.fValScene = fValScene;
        sceneObjects{ss}.multiSceneMeta.stages{stage}.fVal = fVal;
        sceneObjects{ss}.multiSceneMeta.stages{stage}.penalty = penalty;
        sceneObjects{ss}.multiSceneMeta.fVal = fVal;
        sceneObjects{ss}.multiSceneMeta.currentStage = stage;
        
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
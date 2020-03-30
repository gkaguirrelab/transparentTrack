function fVal = multiSceneObjective(x,mySceneObjects,model,multiSceneNorm,verbose)


% Loop over the scenes. Pass the scene params appropriate to that
% scene. Store the objective values.
fValScene = [];
nScenes = length(mySceneObjects);
for ss = 1:nScenes
    fValScene(ss) = mySceneObjects{ss}.updateModel(model.subX(x,ss));
end

% Take the norm of the scene objectives
fVal = norm(fValScene,multiSceneNorm);

% BADS can't handle Inf in the objective, so replace with real max
fVal = min([fVal realmax]);

% Apply the regularization to penalize changes in the camera depth
fVal = fVal * model.penalty(x);

% Each sceneObject has a multiSceneMeta property that is used to stash
% information about the search progress. First determine if we need to
% update the meta information.
updateMetaFlag = false;
if isempty(mySceneObjects{1}.multiSceneMeta)
    updateMetaFlag = true;
else
    if fVal < mySceneObjects{1}.multiSceneMeta.fVal
        updateMetaFlag = true;
    end
end

% We either are initializing or updating the multiSceneMeta information.
if updateMetaFlag
    multiSceneMeta.x = x;
    multiSceneMeta.fVal = fVal;
    multiSceneMeta.fValScene = fValScene;
    multiSceneMeta.penalty = model.penalty(x);
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
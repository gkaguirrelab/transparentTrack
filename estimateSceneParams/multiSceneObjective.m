function fVal = multiSceneObjective(x,mySceneObjects,nEyeParams,nSceneParams,depthChangePenaltyMultiplier,multiSceneNorm,verbose)


% Loop over the scenes. Pass the scene params appropriate to that
% scene. Store the objective values.
fValScene = [];
nScenes = length(mySceneObjects);
for ss = 1:nScenes
    sceneIdx = (ss-1)*nSceneParams+nEyeParams+1;
    subX = x([1:nEyeParams,sceneIdx:sceneIdx+nSceneParams-1]);
    fValScene(ss) = mySceneObjects{ss}.calcError(subX);
end

% Take the norm of the scene objectives
fVal = norm(fValScene,multiSceneNorm);

% BADS can't handle Inf in the objective, so replace with real max
fVal = min([fVal realmax]);

% Apply the regularization to penalize changes in the camera depth
fVal = fVal * depthChangePenaltyMultiplier(x);

% Report the status to the screen
if verbose
    str = 'fVal = %2.2d; x = [ ';
    for pp = 1:length(x)-1
        str = [str '%2.2f, '];
    end
    str = sprintf([str '%2.2f ]\n'],fVal,x);
    fprintf(str);
end

end
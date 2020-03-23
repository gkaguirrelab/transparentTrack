function fVal = multiSceneObjective(x,mySceneObjects,nEyeParams,nSceneParams)

    nScenes = length(mySceneObjects);

    fValScene = [];
    for ss = 1:nScenes
        sceneIdx = (ss-1)*nSceneParams+nEyeParams+1;
        subX = x([1:nEyeParams,sceneIdx:sceneIdx+nSceneParams-1]);
        fValScene(ss) = mySceneObjects{ss}.calcError(subX);
    end

    fVal = norm(fValScene);

end
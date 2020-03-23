function fVal = multiSceneObjective(x,mySceneObjects,nEyeParams,nSceneParams, verbose)

    nScenes = length(mySceneObjects);

    fValScene = [];
    for ss = 1:nScenes
        sceneIdx = (ss-1)*nSceneParams+nEyeParams+1;
        subX = x([1:nEyeParams,sceneIdx:sceneIdx+nSceneParams-1]);
        fValScene(ss) = mySceneObjects{ss}.calcError(subX);
    end

    fVal = norm(fValScene);
    
    if verbose
        str = 'fVal = %2.2d; x = [ ';
        for pp = 1:length(x)-1
            str = [str '%2.2f, '];
        end
        str = sprintf([str '%2.2f ]\n'],fVal,x);
        fprintf(str);
    end

end
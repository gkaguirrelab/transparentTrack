function fVal = calcError(obj, x)

% If the stored is is not empty, and passed x is the same as the current
% state of the model, just return the current fVal
if ~isempty(obj.x)
    if all(obj.x == x)
        fVal = obj.fVal;
        return
    end
end

% Update the sceneGeometry
obj.sceneGeometry = updateSceneGeometry( obj.sceneGeometry, x );

% Calculate the error
[ fVal, modelEyePose, modelPupilEllipse, modelGlintCoord, modelPoseGaze, modelVecGaze, poseRegParams, vectorRegParams, rawErrors] = ...
    calcGlintGazeError( obj.sceneGeometry, obj.args{:}, obj.keyVals{:} );

% Store the parameters and fVal
obj.x = x;
obj.fVal = fVal;

% Store all the other model components
obj.modelEyePose = modelEyePose;
obj.modelPupilEllipse = modelPupilEllipse;
obj.modelGlintCoord = modelGlintCoord;
obj.modelPoseGaze = modelPoseGaze;
obj.modelVecGaze = modelVecGaze;
obj.poseRegParams = poseRegParams;
obj.vectorRegParams = vectorRegParams;
obj.rawErrors = rawErrors;

% Save the fixation 
obj.fixationEyePose = poseRegParams.t;
obj.screenTorsion = poseRegParams.theta;
obj.screenRotMat = poseRegParams.R;


end





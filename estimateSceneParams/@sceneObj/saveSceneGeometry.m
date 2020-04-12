function saveSceneGeometry(obj,fileNameSuffix)
% Saves the current state of the sceneGeometry in the object
%
% Syntax:
%  obj.saveSceneGeometry(fileNameSuffix)
%
% Description:
%   Save the sceneGeometry with the current state of the parameters to
%   disk, using the filepath provided by the videoStemName input to the
%   sceneObj class. The sceneGeometry to be saved is created by specifying
%   a set of key-value arguments that are then passsed to the
%   createSceneGeometry function.
%
% Inputs:
%   fileNameSuffix        - Char vector. A suffix for the saved file name.
%                           The calling function might pass '1', '2', etc
%                           to label different stages of the evolution of
%                           the scene object.
%
% Outputs:
%   none
%

% Get the model parameters and model description
x = obj.x;
model = obj.model;

% The sceneGeometry (and eye) varargin that were passed as part of the
% instantiation of the sceneObj.
sceneGeometryVarargin = obj.setupArgs;

% keys and values to update
keys = {...
    'kvals',...
    'rotationCenterScalers',...
    'primaryPosition', ...
    'cameraTorsion',...
    'cameraTranslation',...
    'poseRegParams',...
    'vecRegParams' ...
    };
values = {...
    x(model.func.fieldSetIdx('eye','kvals')), ...
    x(model.func.fieldSetIdx('eye','rotationCenterScalers')), ...
    x(model.func.fieldSetIdx('scene','primaryPosition')), ...
    x(model.func.fieldParamIdx('scene','torsion')), ...
    x(model.func.fieldSetIdx('scene','translation'))', ...
    obj.poseRegParams, ...
    obj.vecRegParams ...
    };

% Loop through the keys and either update or add the value to the
% sceneGeometry arguments
for kk = 1:length(keys)
    idx = find(strcmp(sceneGeometryVarargin,keys{kk}),1);
    if isempty(idx)
        sceneGeometryVarargin = [sceneGeometryVarargin, keys{kk}, values{kk}];
    else
        sceneGeometryVarargin(idx+1)=values(kk);
    end
end

% Create a new sceneGeometry with the updated key-values
sceneGeometry = createSceneGeometry(sceneGeometryVarargin{:});

% Add the meta data
sceneGeometry.meta.estimateSceneParams.x = obj.x;
sceneGeometry.meta.estimateSceneParams.xHead = obj.x(model.func.fieldSetIdx('head','all'));
sceneGeometry.meta.estimateSceneParams.xEye = obj.x(model.func.fieldSetIdx('eye','all'));
sceneGeometry.meta.estimateSceneParams.xScene = obj.x(model.func.fieldSetIdx('scene','all'));
sceneGeometry.meta.estimateSceneParams.fVal = obj.fVal;
sceneGeometry.meta.estimateSceneParams.sceneGeometryVarargin = sceneGeometryVarargin;
sceneGeometry.meta.estimateSceneParams.obj = obj;
% sceneGeometry.meta.estimateSceneParams.modelEyePose = obj.modelEyePose;
% sceneGeometry.meta.estimateSceneParams.modelPupilEllipse = obj.modelPupilEllipse;
% sceneGeometry.meta.estimateSceneParams.modelGlintCoord = obj.modelGlintCoord;
% sceneGeometry.meta.estimateSceneParams.modelPoseGaze = obj.modelPoseGaze;
% sceneGeometry.meta.estimateSceneParams.modelVecGaze = obj.modelVecGaze;
% sceneGeometry.meta.estimateSceneParams.poseRegParams = obj.poseRegParams;
% sceneGeometry.meta.estimateSceneParams.vecRegParams = obj.vecRegParams;
% sceneGeometry.meta.estimateSceneParams.rawErrors = obj.rawErrors;

% Save the sceneGeometry file
sceneGeometryFileName = [obj.videoStemName '_sceneGeometry' fileNameSuffix '.mat'];
save(sceneGeometryFileName,'sceneGeometry');


end





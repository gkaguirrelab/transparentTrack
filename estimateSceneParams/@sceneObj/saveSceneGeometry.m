function saveSceneGeometry(obj,fileNameSuffix)

% keys and values to update
keys = {...
    'kvals',...
    'rotationCenterScalers',...
    'primaryPosition', ...
    'cameraTorsion',...
    'cameraTranslation',...
    'fixationEyePose',...
    'screenTorsion',...
    'screenRotMat',...
    };
values = {...
    obj.x(1:5), ...
    obj.x(6:7), ...
    obj.x(8:9), ...
    obj.x(10), ...
    obj.x(11:13)', ...
    obj.fixationEyePose, ...
    obj.screenTorsion, ...
    obj.screenRotMat, ...
    };

% Loop through the keys and either update or add
sceneGeometryVarargin = obj.setupArgs;
for kk = 1:length(keys)
    idx = find(strcmp(sceneGeometryVarargin,keys{kk}),1);
    if isempty(idx)
        sceneGeometryVarargin = [sceneGeometryVarargin, keys{kk}, values{kk}];
    else
        sceneGeometryVarargin(idx+1)=values(kk);
    end
end


% Create a new sceneGeometry with the update key-values
sceneGeometry = createSceneGeometry(sceneGeometryVarargin{:});

% Update the meta data
sceneGeometry.meta.estimateSceneParams.p = obj.meta;
sceneGeometry.meta.estimateSceneParams.x = obj.x;
sceneGeometry.meta.estimateSceneParams.fVal = obj.fVal;
sceneGeometry.meta.estimateSceneParams.sceneGeometryVarargin = sceneGeometryVarargin;
sceneGeometry.meta.estimateSceneParams.modelEyePose = obj.modelEyePose;
sceneGeometry.meta.estimateSceneParams.modelPupilEllipse = obj.modelPupilEllipse;
sceneGeometry.meta.estimateSceneParams.modelGlintCoord = obj.modelGlintCoord;
sceneGeometry.meta.estimateSceneParams.modelPoseGaze = obj.modelPoseGaze;
sceneGeometry.meta.estimateSceneParams.modelVecGaze = obj.modelVecGaze;
sceneGeometry.meta.estimateSceneParams.poseRegParams = obj.poseRegParams;
sceneGeometry.meta.estimateSceneParams.vectorRegParams = obj.vectorRegParams;
sceneGeometry.meta.estimateSceneParams.rawErrors = obj.rawErrors;

% Save the sceneGeometry file
sceneGeometryFileName = [obj.videoStemName '_sceneGeometry' fileNameSuffix '.mat'];
save(sceneGeometryFileName,'sceneGeometry');



end





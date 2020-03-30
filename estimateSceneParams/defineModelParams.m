function model = defineModelParams(nScenes)
%% Define model parameters

%% Head
% These parameters adjust the relative camera position vectors derived from
% measurement of head motion during scanning
model.head.labels = {'timeShift','azi','ele','torsion'};
model.head.nParams = length(model.head.labels);
model.head.setLabels = {'phaseAndRotation'};
model.head.setIdx = {1:4};
model.head.idxMap = @(idx) idx;

%% Eye
% These parameters adjust biometric properties of the model eye
model.eye.labels = {'K1','K2','torsion','tilt','tip','joint','diff'};
model.eye.nParams = length(model.eye.labels);
model.eye.setLabels = {'cornea','rotationCenter'};
model.eye.setIdx = {1:5, 6:7};
model.eye.idxMap = @(idx) model.head.nParams+idx;

%% Scene
% These parameters adjust the position of the camera within a scene, which
% also adjusts the primary position of the eye
model.scene.labels = {'pp_azi','pp_ele','torsion','horiz','vert','depth'};
model.scene.nParams = length(model.scene.labels);
model.scene.nScenes = nScenes;
model.scene.setLabels = {'primaryPosition','cameraPosition'};
model.scene.setIdx = {1:2, 3:6};

% An anonymous function that expands the input index vector [a, b, ...]
% into a vector given k eye+relCamPosParams, s sceneParams, and n scenes:
%	[ e+a+(s*0), e+b+(s*0), e+a+(s*1), e+b+(s*1), ... e+a+(s*(n-1)), e+b+(s*(n-1)) ]
model.scene.idxMap = @(idx) model.head.nParams + model.eye.nParams + repmat((0:model.scene.nScenes-1)*model.scene.nParams,1,length(idx)) + ...
    cell2mat(arrayfun(@(x) repmat(x,1,model.scene.nScenes),idx,'UniformOutput',false));

% Arrange the sets into search strages for different search strategies
model.stages.gazeCal = { ...
    {'eye.rotationCenter','scene.cameraPosition'},...
    {'eye.rotationCenter','scene.cameraPosition', 'eye.cornea', 'scene.primaryPosition'} };

model.stages.sceneSync = { ...
    {'scene.cameraPosition'},...
    {'scene.cameraPosition', 'scene.primaryPosition', 'head.phaseAndRotation' } };

% The head and eye parameters are shared by all scenes, and each scene gets
% its own set of scene parameters, yielding this many total parameters in
% the optimization.
model.nParams = model.head.nParams + model.eye.nParams + model.scene.nParams * model.scene.nScenes;


%% subX
% This function returns the indices for the full set of parameters for each
% scene
model.sceneParamStart = @(sceneIdx) (sceneIdx-1)*model.scene.nParams+model.head.nParams+model.eye.nParams+1;
model.subX = @(x,sceneIdx) x([1:(model.head.nParams+model.eye.nParams),model.sceneParamStart(sceneIdx):model.sceneParamStart(sceneIdx)+model.scene.nParams-1]);


%% Depth penalty
% A regularization function that penalizes changes in depth from the x0
% values
cameraDepthTransSet = model.scene.idxMap(find(strcmp('depth',model.scene.labels)));
model.genericPenalty = @(x,x0,w) (1 + w * norm( (x(cameraDepthTransSet) - x0(cameraDepthTransSet)) ./ x0(cameraDepthTransSet) ))^2;


%% Non-linear constraint
% A non-linear constraint for the BADS search that requires first value of
% the corneal curvature (K1) to be less than the second value (K2) Note
% that NONBCON takes a matrix input, which is why we perform this
% calculation over the first dimension.
model.nonbcon = @(x) x(:,model.eye.idxMap(1)) > x(:,model.eye.idxMap(2));


end
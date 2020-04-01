function model = defineModelParams(nScenes, modelIn, verbose)
%% Define model parameters


%% Head
% These parameters adjust the relative camera position vectors derived from
% measurement of head motion during scanning
model.head.x0 = [0, 0, 0, 0];
model.head.bounds = [30, 30, 30, 30];
model.head.paramLabels = {'timeShift','azi','ele','torsion'};
model.head.units = {'seconds','deg','deg','deg'};
model.head.nParams = length(model.head.paramLabels);
model.head.setLabels = {'phaseAndRotation','all'};
model.head.setIdx = {1:4,1:4};
model.head.idxMap = @(idx) idx;
model.head.idxMultiScene = @(idx) idx;


%% Eye
% These parameters adjust biometric properties of the model eye
model.eye.x0 = [44.2410, 45.6302, 0, 2.5000, 0, 1, 1];
model.eye.bounds = [5, 5, 90, 5, 5, 0.25, 0.15];
model.eye.paramLabels = {'K1','K2','torsion','tilt','tip','joint','diff'};
model.eye.units = {'diopters','diopters','deg','deg','deg','proportion','proportion'};
model.eye.nParams = length(model.eye.paramLabels);
model.eye.setLabels = {'kvals','rotationCenterScalers','all'};
model.eye.setIdx = {1:5, 6:7, 1:7};
model.eye.idxMap = @(idx) model.head.nParams+idx;
model.eye.idxMultiScene = @(idx) idx;

%% Scene
% These parameters adjust the position of the camera within a scene, which
% also adjusts the primary position of the eye
model.scene.x0 = [0 0 0 0 0 120];
model.scene.bounds = [10 10 20 20 20 20];
model.scene.paramLabels = {'pp_azi','pp_ele','torsion','horiz','vert','depth'};
model.scene.units = {'deg','deg','deg','mm','mm','mm'};
model.scene.nParams = length(model.scene.paramLabels);
model.scene.nScenes = nScenes;
model.scene.setLabels = {'primaryPosition','cameraPosition', 'translation', 'all'};
model.scene.setIdx = {1:2, 3:6, 4:6, 1:6};
model.scene.idxMap =  @(idx) model.head.nParams+model.eye.nParams+idx;
% An anonymous function that expands the input index vector [a, b, ...]
% into a vector given k eye+head, s sceneParams, and n scenes:
%	[ e+a+(s*0), e+b+(s*0), e+a+(s*1), e+b+(s*1), ... e+a+(s*(n-1)), e+b+(s*(n-1)) ]
% This is used to map a given choice of scene parameters to a multi-scene
% search.
model.scene.idxMultiScene = @(idx) repmat((0:model.scene.nScenes-1)*model.scene.nParams,1,length(idx)) + ...
    cell2mat(arrayfun(@(x) repmat(x,1,model.scene.nScenes),idx,'UniformOutput',false));

%% Stages
% Arrange the sets into search strages for different search strategies
model.stages.gazeCal = { ...
    {'eye.rotationCenterScalers','scene.cameraPosition'},...
    {'eye.rotationCenterScalers','scene.cameraPosition', 'eye.kvals', 'scene.primaryPosition'} };

model.stages.sceneSync = { ...
    {'scene.cameraPosition'},...
    {'scene.cameraPosition', 'scene.primaryPosition', 'head.phaseAndRotation' } };

% The head and eye parameters are shared by all scenes, and each scene gets
% its own set of scene parameters, yielding this many total parameters in
% the optimization.
model.nParams = model.head.nParams + model.eye.nParams + model.scene.nParams * model.scene.nScenes;


%% Substitute passed model inputs for defaults
% We cannot currently overload the functions with the model input as the
% re-ordering of the fields at the time of the merge breaks the function
% linking. Sorry.
if ~isempty(modelIn)
    model = mergestruct(model, modelIn);
end


%% Functions

% Return the indices for a given field (head, eye, scene) and param label
model.func.fieldSetIdx = @(field,setLabel) model.(field).idxMap(model.(field).setIdx{strcmp(model.(field).setLabels,setLabel)});
model.func.fieldParamIdx = @(field,paramLabel) model.(field).idxMap(find(strcmp(model.(field).paramLabels,paramLabel)));


%% subX
% This function returns the indices for the full set of parameters for each
% scene
model.func.sceneParamStart = @(sceneIdx) (sceneIdx-1)*model.scene.nParams+model.head.nParams+model.eye.nParams+1;
model.func.subX = @(x,sceneIdx) x([1:(model.head.nParams+model.eye.nParams),model.func.sceneParamStart(sceneIdx):model.func.sceneParamStart(sceneIdx)+model.scene.nParams-1]);


%% Depth penalty
% A regularization function that penalizes changes in depth from the x0
% values
cameraDepthTransSet = model.scene.idxMultiScene(model.func.fieldParamIdx('scene','depth'));
model.func.genericPenalty = @(x,x0,w) (1 + w * norm( (x(cameraDepthTransSet) - x0(cameraDepthTransSet)) ./ x0(cameraDepthTransSet) ))^2;


%% Non-linear constraint
% A non-linear constraint for the BADS search that requires first value of
% the corneal curvature (K1) to be less than the second value (K2) Note
% that NONBCON takes a matrix input, which is why we perform this
% calculation over the first dimension.
model.func.nonbcon = @(x) x(:,model.func.fieldParamIdx('eye','K1')) > x(:,model.func.fieldParamIdx('eye','K2'));


%% Assemble full x0 and bounds

% If model.scene.x0 is a cell array, then we have been given a different
% set of x0 parameters for each scene.
if iscell(model.scene.x0)
    model.x0 = [model.head.x0, model.eye.x0, cell2mat(model.scene.x0)];
else
    model.x0 = [model.head.x0, model.eye.x0, repmat(model.scene.x0, 1, nScenes)];
end

% And the bounds
if iscell(model.scene.x0)
    model.bounds = [model.head.bounds, model.eye.bounds, cell2mat(model.scene.bounds)];
else
    model.bounds = [model.head.bounds, model.eye.bounds, repmat(model.scene.bounds, 1, nScenes)];
end



end




%% LOCAL FUNCTIONS

function into = mergestruct(into, from)
%MERGESTRUCT merge all the fields of scalar structure from into scalar structure into
validateattributes(from, {'struct'}, {'scalar'});
validateattributes(into, {'struct'}, {'scalar'});
fns = fieldnames(from);
for fn = fns.'
    if isstruct(from.(fn{1})) && isfield(into, fn{1})
        %nested structure where the field already exist, merge again
        into.(fn{1}) = mergestruct(into.(fn{1}), from.(fn{1}));
    else
        %non structure field, or nested structure field that does not already exist, simply copy
        into.(fn{1}) = from.(fn{1});
    end
end
end

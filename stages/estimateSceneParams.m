% x = [7.5817   -1.0905    0.5558  134.0159    0.8400    0.9600    1.1500  0.9903]
% fVal = 1.7198

% housekeeping
clear
close all

% Start the parpool
startParpool(2,true);

% Load the material
load('/Users/aguirre/Dropbox (Aguirre-Brainard Lab)/TOME_processing/session2_spatialStimuli/TOME_3015/032417/EyeTracking/GazeCal03_correctedPerimeter.mat')
load('/Users/aguirre/Dropbox (Aguirre-Brainard Lab)/TOME_processing/session2_spatialStimuli/TOME_3015/032417/EyeTracking/GazeCal03_glint.mat')
load('/Users/aguirre/Dropbox (Aguirre-Brainard Lab)/TOME_processing/session2_spatialStimuli/TOME_3015/032417/EyeTracking/GazeCal03_pupil.mat')
load('/Users/aguirre/Dropbox (Aguirre-Brainard Lab)/TOME_processing/session2_spatialStimuli/TOME_3015/032417/EyeTracking/GazeCal03_sceneGeometry.mat')

% Reduce the material to just the gaze target frames
frames = sceneGeometry.meta.estimateSceneParams.ellipseArrayList;
perimeter.data = perimeter.data(frames);
gazeTargets = sceneGeometry.meta.estimateSceneParams.fixationTargetArray;
ellipseRMSE = pupilData.initial.ellipses.RMSE(frames);
glintData.X = glintData.X(frames); glintData.Y = glintData.Y(frames);

% Define BADS search options
options = bads('defaults');          % Get a default OPTIONS struct
options.Display = 'off';             % Silence display output
options.UncertaintyHandling = 0;     % The objective is deterministic

% Define the complete objective function in which we derive eyePose
myObjComplete = @(x) calcGlintGazeError( updateSceneGeometry( sceneGeometry, x ), perimeter, gazeTargets, ellipseRMSE, glintData, [] );

% Define an add plots function
nStages = 4;
addPlotsWrap = @(idx,x) addSubPlots(idx,x,nStages,sceneGeometry,perimeter, gazeTargets, ellipseRMSE, glintData);

% Start a figure
figure


%% STAGE 1 -- TORSION / TRANSLATION SEARCH
% Perform an initial, iterated search, locking parameters for camera
% distance, eye rotation, and corneal curvature.

x = [0, 0, 0, 130, 1, 1, 1, 1];
bound = [20, 10, 10, 0, 0, 0, 0, 0];
lb = x - bound;
ub = x + bound;
lbp = x - bound./2;
ubp = x + bound./2;

xLast = x;
fValLast = realmax;
stillSearching  = true;
fprintf('Stage 1...');
while stillSearching
    % obtain the modelEyePose
    [ ~, modelEyePose] = calcGlintGazeError( updateSceneGeometry( sceneGeometry, x ), perimeter, gazeTargets, ellipseRMSE, glintData, []);
    % Update the objective
    myObj = @(x) calcGlintGazeError( updateSceneGeometry( sceneGeometry, x ), perimeter, gazeTargets, ellipseRMSE, glintData, modelEyePose );
    % Perform an initial fit to optimize glint
    x = bads(myObj,x,lb,ub,lbp,ubp,[],options);
    % The objective we care about is for the complete objective in which we
    % calculate the eye pose for this x.
    fVal = myObjComplete(x);
    % Update progress
    fprintf([num2str(fVal) '...']);
    % Check for done searching
    if fVal >= fValLast
        x = xLast;
        fVal = fValLast;
        stillSearching = false;
        fprintf([num2str(fVal) '\n']);
        x
    else
        xLast = x;
        fValLast = fVal;
    end
end

% Update the plot
addPlotsWrap(1,x);


%% STAGE 2 -- ROTATION CENTER SEARCH
% Search over the eye rotation center

% Hard and plausible bounds
lb = [x(1:4), 0.75, 0.75, x(7:8)];
ub = [x(1:4), 1.25, 1.25, x(7:8)];
lbp = [x(1:4), 0.85, 0.85, x(7:8)];
ubp = [x(1:4), 1.15, 1.15, x(7:8)];

% The objective considers only gazeError
myObj = @(x) calcGlintGazeError( updateSceneGeometry( sceneGeometry, x ), perimeter, gazeTargets, ellipseRMSE, glintData, [] );

% Search
fprintf('Stage 2...');
[x, fVal] = bads(myObj,x,lb,ub,lbp,ubp,[],options);
fprintf([num2str(fVal) '\n']);
x

% Update the plot
addPlotsWrap(2,x);


%% STAGE 3 -- TRANSLATION AND CURVATURE SEARCH
% Lock the rotation centers, search over translation and corneal curvature
bound = [abs(x(1:3).*0.25), 0, 0, 0, x(7:8).*0.25];
lb = x - bound;
ub = x + bound;
lbp = x - bound./2;
ubp = x + bound./2;

xLast = x;
fValLast = realmax;
stillSearching  = true;
fprintf('Stage 3...');
while stillSearching
    % obtain the modelEyePose
    [ ~, modelEyePose] = calcGlintGazeError( updateSceneGeometry( sceneGeometry, x ), perimeter, gazeTargets, ellipseRMSE, glintData, []);
    % Update the objective
    myObj = @(x) calcGlintGazeError( updateSceneGeometry( sceneGeometry, x ), perimeter, gazeTargets, ellipseRMSE, glintData, modelEyePose );
    % Perform the search
    x = bads(myObj,x,lb,ub,lbp,ubp,[],options);
    % The objective we care about is for the complete objective
    fVal = myObjComplete(x);
    % Update progress
    fprintf([num2str(fVal) '...']);
    if fVal >= fValLast
        x = xLast;
        fVal = fValLast;
        stillSearching = false;
        fprintf([num2str(fVal) '\n']);
        x
    else
        xLast = x;
        fValLast = fVal;
    end
end

% Update the plot
addPlotsWrap(3,x);



%% STAGE 4 -- COMPLETE SEARCH

% Hard and plausible bounds
lb  = x./(0.90.^-sign(x));
lbp = x./(0.95.^-sign(x));
ubp = x./(1.05.^-sign(x));
ub  = x./(1.10.^-sign(x));

% Standard objective
myObj = @(x) calcGlintGazeError( updateSceneGeometry( sceneGeometry, x ), perimeter, gazeTargets, ellipseRMSE, glintData, [] );

% Search
fprintf('Stage 4...');
[x, fVal] = bads(myObj,x,lb,ub,lbp,ubp,[],options);
fprintf([num2str(fVal) '\n']);
x

% Update the plot
addPlotsWrap(4,x);


%% Update the sceneGeometry 

% The fitted sceneGeometry
f = updateSceneGeometry( sceneGeometry, x );

% The varargin originally used to create the sceneGeometry
sceneGeometryVarargin = sceneGeometry.meta.createSceneGeometry.varargin;

% keys and values to update
keys = {...
    'cameraTorsion',...
    'cameraTranslation',...
    'rotationCenters',...
    'measuredCornealCurvature'};
values = {...
    f.cameraPosition.torsion, ...
    f.cameraPosition.translation, ...
    f.eye.rotationCenters, ...
    f.eye.meta.measuredCornealCurvature};

% Loop through the keys and either update or add
for kk = 1:length(keys)
    idx = find(strcmp(sceneGeometryVarargin,keys{kk}),1);
    if isempty(idx)
        sceneGeometryVarargin = [sceneGeometryVarargin, keys{kk}, values{kk}];
    else
        sceneGeometryVarargin(idx+1)=values(kk);
    end
end

% Save the metaData of the original sceneGeometry creation
metaCreate = sceneGeometry.meta.createSceneGeometry;

% Create a new sceneGeometry with the update key-values
newSceneGeometry = createSceneGeometry(sceneGeometryVarargin{:});

% Update and move the meta data around
newSceneGeometry.meta.estimateSceneParams.create = newSceneGeometry.meta.createSceneGeometry;
newSceneGeometry.meta.createSceneGeometry = metaCreate;
newSceneGeometry.meta.estimateSceneParams.x = x;
    

%%%%%%%%%%%% LOCAL FUNCTIONS

function addSubPlots(idx,x,nStages,sceneGeometry,perimeter,gazeTargets, ellipseRMSE, glintData)

[ ~, ~, modelGlint, modelPoseGaze, modelVecGaze] = ...
    calcGlintGazeError( updateSceneGeometry( sceneGeometry, x ), perimeter, gazeTargets, ellipseRMSE, glintData, []);

nCols = 3;

subplot(nStages,nCols,(idx-1)*nCols+1)
plot(gazeTargets(1,:),gazeTargets(2,:),'ok'); hold on;
plot(modelPoseGaze(1,:),modelPoseGaze(2,:),'xr'); hold on;
ylim([-10 10])
axis equal
title(['Stage ' num2str(idx) ' -- poseGaze']);


subplot(nStages,nCols,(idx-1)*nCols+2)
plot(gazeTargets(1,:),gazeTargets(2,:),'ok'); hold on;
plot(modelVecGaze(1,:),modelVecGaze(2,:),'xr'); hold on;
ylim([-10 10])
axis equal
title(['Stage ' num2str(idx) ' -- vecGaze']);


subplot(nStages,nCols,(idx-1)*nCols+3)
plot(glintData.X,glintData.Y,'ok'); hold on;
plot(modelGlint.X,modelGlint.Y,'xr'); hold on;
axis equal
title(['Stage ' num2str(idx) ' -- glint']);

drawnow

end

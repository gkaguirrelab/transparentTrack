function sceneGeometry = estimateSceneParams(pupilFileName, perimeterFileName, glintFileName, sceneGeometryFileName, varargin)
% Estimate scene geometry parameters from pupil and glint measurements
%
% Syntax:
%  sceneGeometry = estimateSceneParams(pupilFileName, perimeterFileName, glintFileName, sceneGeometryFileName, varargin)
%
% Description:
%   The appearance of the eye in a camera image is influenced by several
%   parameters, including biometric properties of the eye and the position
%   of the camera. This function uses a set of observations of the pupil
%   and the glint to estimate these scene parameters. The accuracy of the
%   estimate is much improved by supplying the location of visual targets
%   that the eye was fixated upon during each of the observations.
%
% Inputs:
%	pupilFileName
%   perimeterFileName
%   glintFileName         - Full path to these files.
%   sceneGeometryFileName - Full path to the file in which the
%                           sceneGeometry data should be saved
%
% Optional key/value pairs (display and I/O):
%  'verbose'              - Logical. Default false.
%  'pupilFileToVideoSuffixSwitch' - Cell array that provides the suffix
%                           of the pupilData file and the suffix of the
%                           corresponding fit video file. This way, the fit
%                           video corresponding to the passed pupilData
%                           file can be found and used to create the
%                           ellipse array montage plot.
%
% Optional key/value pairs (flow control)
%  'useParallel'          - If set to true, use the MATLAB parallel pool
%  'nWorkers'             - Specify the number of workers in the parallel
%                           pool. If undefined the default number will be
%                           used.
%
% Optional key/value pairs (environment)
%  'tbSnapshot'           - This should contain the output of the
%                           tbDeploymentSnapshot performed upon the result
%                           of the tbUse command. This documents the state
%                           of the system at the time of analysis.
%  'timestamp'            - AUTOMATIC; The current time and date
%  'username'             - AUTOMATIC; The user
%  'hostname'             - AUTOMATIC; The host
%
% Optional key/value pairs (analysis)
%  'sceneParamsLB/UB'     - 6x1 vector. Hard upper and lower bounds. Should
%                           reflect the physical limits of the measurement.
%  'sceneParamsLBp/UBp'   - 6x1 vector. Plausible upper and lower bounds.
%                           Where you think the translation vector solution
%                           is likely to be.
%  'eyePoseLB/UB'         - 1x4 vector. Upper / lower bounds on the eyePose
%                           [azimuth, elevation, torsion, pupil radius].
%                           The torsion value is unusued and is bounded to
%                           zero. Biological limits in eye rotation and
%                           pupil size would suggest boundaries of [±35,
%                           ±25, 0, 0.25-4]. Note, however, that these
%                           angles are relative to the center of
%                           projection, not the primary position of the
%                           eye. Therefore, in circumstances in which the
%                           camera is viewing the eye from an off-center
%                           angle, the bounds will need to be shifted
%                           accordingly.
%  'frameSet'             - A vector of m frame numbers (indexed from 1)
%                           which identify the ellipses to be used for the
%                           estimation of scene geometry. If left empty,
%                           a list of ellipses will be generated.
%  'gazeTargets'          - A 2xm matrix that provides the positions, in
%                           nominal degrees of visual angle, of fixation
%                           targets that correspond to each of the frames
%                           identified in the ellipseArrayList. If defined,
%                           the routine will find the scene geometry that
%                           best aligns the recovered eye poses with the
%                           fixation targets, subject to a translation and
%                           rotation matrix. If left empty, the search will
%                           minimize error in the joint specification of
%                           ellipse centers and shape. If needed, the
%                           visual angle of the stimuli will be adjusted
%                           for min/magnification produced by spectacle
%                           lenses worn by the subject.
%  'fixSpectacleLens'     - Scalar. This parameter is used to handle an
%                           unusual circumstance in which the eye viewing a
%                           fixation array is behind a spectacle lens, but
%                           the eye being modeled for pupil appearance is
%                           not. Setting this parameter causes the routine
%                           to calculate a magnification factor for the
%                           fixation target array, but does not apply this
%                           spectacle to the eye being modeled.
%  'nBADSsearches'        - Scalar. We perform the search over scene params
%                           from a randomly selected starting point within
%                           the plausible bounds. This parameter sets how
%                           many random starting points to try; the best
%                           result is retained. Each search is run on a
%                           separate worker if the parpool is available. If
%                           a value of zero is passed, then a sceneGeometry
%                           file and diagnostic plots are created using the
%                           midpoint of the passed bounds.

%
% Outputs
%	sceneGeometry         - A structure that contains the components of the
%                           projection model.
%
% Examples:
%{
    perimeterFileName = '/Users/aguirre/Dropbox (Aguirre-Brainard Lab)/TOME_processing/session2_spatialStimuli/TOME_3015/032417/EyeTracking/GazeCal03_correctedPerimeter.mat';
    glintFileName = '/Users/aguirre/Dropbox (Aguirre-Brainard Lab)/TOME_processing/session2_spatialStimuli/TOME_3015/032417/EyeTracking/GazeCal03_glint.mat';
    pupilFileName = '/Users/aguirre/Dropbox (Aguirre-Brainard Lab)/TOME_processing/session2_spatialStimuli/TOME_3015/032417/EyeTracking/GazeCal03_pupil.mat';

    load('/Users/aguirre/Dropbox (Aguirre-Brainard Lab)/TOME_processing/session2_spatialStimuli/TOME_3015/032417/EyeTracking/GazeCal03_sceneGeometry.mat')
    frameSet = sceneGeometry.meta.estimateSceneParams.ellipseArrayList;
    gazeTargets = sceneGeometry.meta.estimateSceneParams.fixationTargetArray;
    sceneGeometryFileName = '~/Desktop/foo_sceneGeometry.mat';

    estimateSceneParams(pupilFileName, perimeterFileName, glintFileName, sceneGeometryFileName, 'frameSet', frameSet, 'gazeTargets', gazeTargets)
%}


%% input parser
p = inputParser; p.KeepUnmatched = true;

% Required
p.addRequired('pupilFileName',@ischar);
p.addRequired('perimeterFileName',@ischar);
p.addRequired('glintFileName',@ischar);
p.addRequired('sceneGeometryFileName',@ischar);

% Optional display and I/O params
p.addParameter('verbose',false,@islogical);
p.addParameter('grayVideoName','',@(x)(isempty(x) | ischar(x)));
p.addParameter('pupilFileToVideoSuffixSwitch',{'_pupil.mat','_gray.avi'},@iscell);

% Optional flow control params
p.addParameter('useParallel',false,@islogical);
p.addParameter('nWorkers',[],@(x)(isempty(x) || isnumeric(x)));

% Optional environment parameters
p.addParameter('tbSnapshot',[],@(x)(isempty(x) || isstruct(x)));
p.addParameter('timestamp',char(datetime('now')),@ischar);
p.addParameter('username',char(java.lang.System.getProperty('user.name')),@ischar);
p.addParameter('hostname',char(java.net.InetAddress.getLocalHost.getHostName),@ischar);

% Optional analysis params
p.addParameter('sceneParamsX0',[],@isnumeric);
p.addParameter('eyePoseLB',[-89,-89,0,0.1],@isnumeric);
p.addParameter('eyePoseUB',[89,89,0,4],@isnumeric);
p.addParameter('frameSet',[],@(x)(isempty(x) | isnumeric(x)));
p.addParameter('gazeTargets',[],@(x)(isempty(x) | isnumeric(x)));
p.addParameter('fixSpectacleLens',[],@(x)(isempty(x) | isnumeric(x)));
p.addParameter('nBinsPerDimension',4,@isnumeric);
p.addParameter('badFrameErrorThreshold',2, @isnumeric);

% parse
p.parse(pupilFileName, perimeterFileName, glintFileName, sceneGeometryFileName, varargin{:})




%% Announce we are starting
if p.Results.verbose
    tic
    fprintf(['Estimating scene parameters. Started ' char(datetime('now')) '\n']);
end


%% Create initial sceneGeometry structure
sceneGeometry = createSceneGeometry(varargin{:});


%% Define the fixationTargetArray
gazeTargets = p.Results.gazeTargets;


%% Handle spectacle magnification
% A spectacle lens has the property of magnifying / minifying the visual
% world from the perspective of the eye. This alteration scales the
% apparent visual field positions of the targets and results in a
% concomittant change in eye movement amplitude. Note that while a contact
% lens also has a magnification effect (albeit smaller), the lens rotates
% with the eye. Thus, eye movement amplitude is not altered.
if ~isempty(gazeTargets)
    % Default to no change
    magnification = 1;
    % If fixSpectacleLens is set, use this value to calculate a
    % magnification and apply it
    if ~isempty(p.Results.fixSpectacleLens)
        modVarargin = varargin;
        idx = find(strcmp(modVarargin,'spectacleLens'));
        if ~isempty(idx)
            modVarargin{idx+1} = p.Results.fixSpectacleLens;
        else
            modVarargin = [modVarargin 'spectacleLens' p.Results.fixSpectacleLens];
        end        
        tmpSceneGeometry = createSceneGeometry(modVarargin{:});
        magnification = tmpSceneGeometry.refraction.retinaToCamera.magnification.spectacle;
    else
        % check if there is a spectacle magnification field
        if isfield(sceneGeometry.refraction.retinaToCamera.magnification,'spectacle')
            magnification = sceneGeometry.refraction.retinaToCamera.magnification.spectacle;
        end
    end
    gazeTargets = gazeTargets .* magnification;
end


%% Load the materials
load(pupilFileName,'pupilData');
load(perimeterFileName,'perimeter');
load(glintFileName,'glintData');


%% Restrict the materials to the frameSet

% Handle the frameset
frameSet = p.Results.frameSet;
if isempty(frameSet)
    % Call here to selectFrameSet    
end

% Extract the frames we want
perimeter.data = perimeter.data(frameSet);
ellipseRMSE = pupilData.initial.ellipses.RMSE(frameSet);
glintData.X = glintData.X(frameSet); glintData.Y = glintData.Y(frameSet);

% Assemble these components into the args variable
args = {perimeter, gazeTargets, ellipseRMSE, glintData};


%% Set up the parallel pool
if p.Results.useParallel
    nWorkers = startParpool( p.Results.nWorkers, p.Results.verbose );
else
    nWorkers=0;
end


%% Set up a figure
figure
nStages = 4;
addPlotsWrap = @(idx,x) addSubPlots(idx,x,nStages,sceneGeometry,args{:});


%% Define BADS search options
options = bads('defaults');          % Get a default OPTIONS struct
options.Display = 'off';             % Silence display output
options.UncertaintyHandling = 0;     % The objective is deterministic


%% STAGE 1 -- TORSION / TRANSLATION SEARCH
% Perform an initial, iterated search, locking parameters for camera
% distance, eye rotation, and corneal curvature.

x = [0, 0, 0, 130, 1, 1, 1, 1];

% Set bounds
bound = [20, 10, 10, 0, 0, 0, 0, 0];
lb = x - bound;
ub = x + bound;
lbp = x - bound./2;
ubp = x + bound./2;

% Search
x = iterativeSearch(x,sceneGeometry,args,lb,ub,lbp,ubp, options);
addPlotsWrap(1,x);


%% STAGE 2 -- ROTATION CENTER SEARCH
% Search over the eye rotation center

% Set bounds
lb = [x(1:4), 0.75, 0.75, x(7:8)];
ub = [x(1:4), 1.25, 1.25, x(7:8)];
lbp = [x(1:4), 0.85, 0.85, x(7:8)];
ubp = [x(1:4), 1.15, 1.15, x(7:8)];

% Set the objective
myObj = @(x) calcGlintGazeError( updateSceneGeometry( sceneGeometry, x ), perimeter, gazeTargets, ellipseRMSE, glintData, [] );

% Search
%x = bads(myObj,x,lb,ub,lbp,ubp,[],options);
addPlotsWrap(2,x);


%% STAGE 3 -- TRANSLATION AND CURVATURE SEARCH
% Lock the rotation centers, search over translation and corneal curvature

% Set bounds
bound = [abs(x(1:3).*0.25), 0, 0, 0, x(7:8).*0.25];
lb = x - bound;
ub = x + bound;
lbp = x - bound./2;
ubp = x + bound./2;

% Search
%x = iterativeSearch(x,sceneGeometry,args,lb,ub,lbp,ubp, options);
addPlotsWrap(3,x);


%% STAGE 4 -- COMPLETE SEARCH

% Set bounds
lb  = x./(0.90.^-sign(x));
lbp = x./(0.95.^-sign(x));
ubp = x./(1.05.^-sign(x));
ub  = x./(1.10.^-sign(x));

% Set the objective
myObj = @(x) calcGlintGazeError( updateSceneGeometry( sceneGeometry, x ), perimeter, gazeTargets, ellipseRMSE, glintData, [] );

% Search
%[x, fVal] = bads(myObj,x,lb,ub,lbp,ubp,[],options);
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
    f.eye.cornea.kvals};

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
    
end






%%%%%%%%%%%% LOCAL FUNCTIONS




function [x, fVal] = iterativeSearch(x,sceneGeometry,args,lb,ub,lbp,ubp, options)
% Implements an iterative search for scene parameters
%
% Syntax:
%  [x, fVal] = iterativeSearch(x,sceneGeometry,args,lb,ub,lbp,ubp)
%
% Description:
%   
%   
xLast = x;
fValLast = realmax;
stillSearching  = true;

while stillSearching
    
    % obtain the modelEyePose
    [ ~, modelEyePose] = calcGlintGazeError( updateSceneGeometry( sceneGeometry, x ), args{:}, []);
    % Update the objective
    myObj = @(x) calcGlintGazeError( updateSceneGeometry( sceneGeometry, x ), args{:}, modelEyePose );
    % Perform the search
    x = bads(myObj,x,lb,ub,lbp,ubp,[],options);
    % The objective we care about is for the complete objective
    fVal = calcGlintGazeError( updateSceneGeometry( sceneGeometry, x ), args{:}, [] );
    % Evalaute the results
    if fVal >= fValLast
        x = xLast;
        fVal = fValLast;
        stillSearching = false;
    else
        xLast = x;
        fValLast = fVal;
    end
end

end





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

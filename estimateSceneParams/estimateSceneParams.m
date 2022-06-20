function sceneObjects = estimateSceneParams(videoStemName, frameSet, gazeTargets, varargin)
% Estimate scene geometry and head motion parameters for one or more scenes
%
% Syntax:
%  estimateSceneParams(videoStemName, frameSet, gazeTargets)
%
% Description:
%   The appearance of the eye in a camera image is influenced by several
%   parameters, including biometric properties of the eye and the position
%   of the camera. The relative position of the camera can further be
%   altered by motion of the head over time. This function uses a set of
%   observations of the pupil and the glint to estimate these scene
%   parameters. This may include the location of visual targets that the
%   eye was fixated upon during each of the observations. A set of
%   measurements made with the camera in a fixed location constitutes a
%   "scene". The routine may be supplied with information from more than
%   one scene. If available, a measurement of the relative position of the
%   head over time will be integrated into the prediction.
%
%   The behavior of the search is principally controlled by using
%   'searchStrategy' key. The values are:
%
%   gazeCal   - Expects one or more gaze cal videos, for which a gazeTarget
%               in degrees of visual angle is specified for each frame of
%               observation of the eye. If multiple scenes are specified,
%               the search forces the eye biometric parameters to be the
%               same in all scenes, but allows the scene parameters (camera
%               position) to vary between scenes. The search does not
%               attempt to adjust the effect of head motion upon relative
%               camera position.
%
%   sceneSync - Expects a single video. Typically, one would supply the eye
%               biometric parameters (estimated elsewhere or measured) by
%               passing model.eye.x0. The search is performed across camera
%               position, including adjustments of the effects of head
%               motion upon relative camera position.
%
%   The accuracy of the solution is improved by setting the cameraDepth and
%   cameraTorsion key-values. The cameraDepth is the distance of the camera
%   from the eye (in mm). It is difficult for the solution to distinguish
%   between the depth position of the camera and properties of eye biometry
%   (rotation depth and corneal curvature). The utility estimateCameraDepth
%   may be used to supply this value by reference to the width of the iris.
%   The cameraTorsion is the rotation of the camera with respect to the
%   azimuthal plane of eye rotation. It is difficult for the routine to
%   distinguish between camera rotation and rotation of the cornea about
%   the optical axis of the eye. The utility estimateCameraTorsion may be
%   used to supply this value by reference to the angle between the medial
%   and lateral canthii of the eye. The search includes a regularization
%   that penalizes the solution for departing from the x0 values for camera
%   torsion and depth.
%
% Inputs:
%	videoStemName         - Char vector, or cell array of n char vectors.
%                           Full path to the n video files from which the
%                           scene observations have been derived. The stem
%                           name should omit the "_gray.avi" suffix that is
%                           usually present in the names of these video
%                           files. The routine assumes that files to be
%                           used in the analysis have the suffixes:
%                               {'_correctedPerimeter.mat','_glint.mat',...
%                                '_pupil.mat'}
%                           and optionally {'_relativeCameraPosition.mat'}.
%                           The sceneGeometry output files will be saved at
%                           the same path location with the suffix:
%                               '_sceneGeometry.mat'
%   frameSet              - A 1xm vector, or a cell array of n such vectors
%                           (m may differ in each vector). Each vector
%                           specifies the m frame indices (indexed from 1)
%                           which identify the set of frames from the scene
%                           to guide the search.
%   gazeTargets           - A 2xm matrix, or cell array of n such matrices,
%                           (m may differ in each matrix). Each matrix
%                           provides the positions, in degrees of visual
%                           angle, of fixation targets that correspond to
%                           each of the frames. The visual angle of the
%                           stimuli should be adjusted for minification /
%                           magnification produced by spectacle lenses worn
%                           by the subject prior to being passed to this
%                           routine.
%
% Optional key/value pairs (display and I/O):
%  'verbose'              - Logical. Default false.
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
%  'outputFileSuffix'     - Char vector. A string that is appended to the
%                           file name of the saved sceneGeometry and
%                           diagnostic plots.
%  'searchStrategy'       - Char vector. Controls the search behavior. See
%                           the Description section above.
%  'cameraTorsion'        - Scalar or vector of length n. Rotation of the
%                           camera with respect to the azimuthal plane of
%                           rotation of the eye, in degrees. Used to set
%                           the scene x0.
%  'cameraDepth'          - Scalar or vector of length n. Estimate of the
%                           distance of the nodal point of the camera from
%                           the corneal apex of the eye in mm. Used to
%                           inform the x0 values for the search.
%  'corneaTorsion'        - Scalar. The angle of astigmatism for the cornea
%                           that is used to set the x0 value, perhaps
%                           obtained from keratometry measurement for the
%                           eye to be modeled. This is particularly useful
%                           for angles close to 90 degrees, for which the
%                           model has trouble reaching by search on its
%                           own.
%  'model'                  Structure. Set fields of this structure to
%                           replace the default set of params for the
%                           search. See defineModelParams.m for details.
%  'eyeArgs'              - Cell array. These are key-value pairs to be
%                           used when generating the eye model for this
%                           subject. An example input is:
%                               {'axialLength',23.45,'sphericalAmetropia',-0.5}
%                           Note that these key values:
%                               {'kvals','rotationCenterScalers','primaryPosition'}
%                           should not be passed, as they will be
%                           supplanted by parameters in the search. If you
%                           wish to provide and lock these parameters in
%                           the search, do so by providing the values in
%                           the x0 arrays by passing the model key-value.
%  'sceneArgs'            - Cell array of key-values, or an array of n such
%                           cell arrays. Each array contians key-value
%                           pairs to be used when generating the scene
%                           model for a particular video. An example input
%                           (for 2 scenes ) is:
%                               { {'spectacleLens',-3} ,{'contactLens',-3} }
%  'errorArgs'            - Cell array. These are key-value pairs that are
%                           passed to method 'updateError' of sceneObj.
%                           Typically, one might pass pre-computed values
%                           of poseRegParams, or adjust the eyePose bounds.
%                           Note that the errorReg key-value should be
%                           defined as part of a search strategy within
%                           defineModelParams.m, not passed here.
%
% Outputs
%   sceneObjects          - Cell array of handles to objects of the type
%                           sceneObj.
%
% Examples:
%{
    % Get the DropBox base directory
    dropboxBaseDir = getpref('eyeTrackTOMEAnalysis','dropboxBaseDir');

    % A set of five gazeCal acquisitions, and the set of gazeTargets and
    % corresponding frames
    videoStemName = {fullfile(dropboxBaseDir,'TOME_processing/session2_spatialStimuli/TOME_3015/032417/EyeTracking/GazeCal01'), ...
                    fullfile(dropboxBaseDir,'TOME_processing/session2_spatialStimuli/TOME_3015/032417/EyeTracking/GazeCal02'), ...
                    fullfile(dropboxBaseDir,'TOME_processing/session2_spatialStimuli/TOME_3015/032417/EyeTracking/GazeCal03'), ...
                    fullfile(dropboxBaseDir,'TOME_processing/session2_spatialStimuli/TOME_3015/032417/EyeTracking/GazeCal04'), ...
                    fullfile(dropboxBaseDir,'TOME_processing/session1_restAndStructure/TOME_3015/030117/EyeTracking/GazeCal')};
    gazeTargets = {...
        [ -7, 7, 0, 0, 7, -7, 7, -7 ; 7, 0, 0, 7, 7, 0, -7, -7], ...
    	[ 0, 0, -7, 7, 0, 7, -7, 7, -7 ; -7, 7, 7, -7, 0, 7, -7, 0, 0], ...
        [ -7, 0, -7, 7, 7, 0, -7, 0, 7 ; 0, -7, -7, 0, -7, 7, 7, 0, 7], ...
        [ 7, -7, -7, 0, 0, 7, -7, 7, 0 ; 7, 0, 7, 7, -7, -7, -7, 0, 0], ...
        [ 0, -7, -7, -7, 7, 0, 7, 0, 7 ; -7, 7, -7, 0, -7, 0, 7, 7, 0]};
    frameSet = {...
        [ 697, 1166, 1300, 1431, 1569, 1711, 1910, 2017 ], ...
        [ 621, 818, 1117, 1277, 1487, 1668, 1864, 2004, 2134 ], ...
        [ 679, 884, 1180, 1250, 1571, 1663, 1809, 2004, 2075 ], ...
        [ 573, 710, 921, 1151, 1285, 1483, 1617, 1800, 1878 ], ...
        [ 660, 842, 904, 1077, 1230, 1328, 1477, 1614, 1670 ]};

    % Fixed biometric properties of the eye. We may know something about
    % corneal curvature, but we won't set that here as we wish to estimate
    % those params from the search. If we wanted, we could pass our
    % measurement as an x0 value using the model key-value.
    eyeArgs = {'axialLength',23.45,'sphericalAmetropia',-0.5};

    % Perform the search
    estimateSceneParams(videoStemName, frameSet, gazeTargets, ...
        'searchStrategy','gazeCal','cameraDepth',140,'eyeArgs',eyeArgs);
%}
%{
    % This is a sceneGeometry file created from the analysis of gazeCal
    % acquisitions. We want to use the eye parameters from that solution,
    % including the registration of eye pose to fixation angle.

    % Get the DropBox base directory
    dropboxBaseDir = getpref('eyeTrackTOMEAnalysis','dropboxBaseDir');

    % Get eye biometric information from the source sceneGeometry file
    sourceDir = fullfile(dropboxBaseDir,'TOME_processing/session2_spatialStimuli/TOME_3015/032417/EyeTracking/');
    sourceSceneName = fullfile(sourceDir,'GazeCal01_sceneGeometry.mat');
    load(sourceSceneName,'sceneGeometry');
    model.eye.x0 = sceneGeometry.meta.estimateSceneParams.x(5:11);
    eyeArgs = sceneGeometry.meta.estimateSceneParams.p.eyeArgs;
    errorArgs = { ...
        'poseRegParams',sceneGeometry.meta.estimateSceneParams.poseRegParams,...
        'vectorRegParams',sceneGeometry.meta.estimateSceneParams.vectorRegParams};

    % Get the cameraDepth from the source sceneGeometry file
    cameraDepth = sceneGeometry.meta.estimateSceneParams.x(end);

    % This is the video for which we wish to create a sceneGeometry file.
    % We select frames to guide the search, using a fixation period before
    % the scan and a distributed set of gaze positions after the scan start
    videoStemName = fullfile(sourceDir,'tfMRI_MOVIE_PA_run03');
    [frameSet1, gazeTargets1] = selectFrames.gazePre(videoStemName);
    [frameSet2, gazeTargets2] = selectFrames.gridSpace(videoStemName);
    frameSet = [frameSet1 frameSet2];
    gazeTargets = [gazeTargets1 gazeTargets2];

    % Perform the search
    estimateSceneParams(videoStemName, frameSet, gazeTargets, ...
        'searchStrategy','sceneSync','cameraDepth',cameraDepth,'model',model,...
        'eyeArgs',eyeArgs,'errorArgs',errorArgs);
%}



%% input parser
p = inputParser; p.KeepUnmatched = true; p.PartialMatching = false;

% Required
p.addRequired('videoStemName',@(x)(iscell(x) || ischar(x)));
p.addRequired('frameSet',@(x)(iscell(x) || isnumeric(x)));
p.addRequired('gazeTargets',@(x)(iscell(x) || isnumeric(x)));

% Optional display and I/O params
p.addParameter('verbose',true,@islogical);
p.addParameter('savePlots',true,@islogical);
p.addParameter('saveFiles',true,@islogical);

% Optional flow control params
p.addParameter('useParallel',true,@islogical);
p.addParameter('nWorkers',[],@(x)(isempty(x) || isnumeric(x)));

% Optional environment parameters
p.addParameter('tbSnapshot',[],@(x)(isempty(x) || isstruct(x)));
p.addParameter('timestamp',char(datetime('now')),@ischar);
p.addParameter('username',char(java.lang.System.getProperty('user.name')),@ischar);
p.addParameter('hostname',char(java.net.InetAddress.getLocalHost.getHostName),@ischar);

% Optional analysis params
p.addParameter('outputFileSuffix','',@ischar);
p.addParameter('searchStrategy','gazeCal',@ischar);
p.addParameter('cameraTorsion',0,@(x)(isscalar(x) || isvector(x)));
p.addParameter('cameraDepth',120,@(x)(isscalar(x) || isvector(x)));
p.addParameter('corneaTorsion',0,@isscalar);
p.addParameter('model',[],@isstruct);
p.addParameter('eyeArgs',{},@iscell);
p.addParameter('sceneArgs',{},@iscell);
p.addParameter('errorArgs',{},@iscell);
p.addParameter('glintData',{},@iscell);
p.addParameter('perimeter',{},@iscell);
p.addParameter('relativeCameraPosition',{},@iscell);

% Optional data input


% parse
p.parse(videoStemName, frameSet, gazeTargets, varargin{:})


%% Check inputs

% Catch some possible input errors
if iscell(videoStemName)
    % If the videoNameStem is a cell array, make sure that all inputs are cell
    % arrays
    if ~iscell(frameSet) || ~iscell(gazeTargets)
        error('estimateSceneParams:inputTypeError','If one required input is a cell array, they must all be cell arrays')
    end
    if ~isequal(length(videoStemName),length(frameSet),length(gazeTargets))
        error('estimateSceneParams:inputTypeError','The cell array inputs must be of equal length')
    end
else
    % The other valid type for the videoNameStem is a single char vector
    if ~ischar(videoStemName)
        error('estimateSceneParams:inputTypeError','videoStemName must be either a char vector or a cell array')
    end
    % In this case, the frameSet and gazeTargets must be numeric, with the same
    % size first dimension
    if ~isnumeric(frameSet) || ~isnumeric(gazeTargets)
        error('estimateSceneParams:inputTypeError','If the videoStemName is a char vector, the other required inputs must be numeric')
    end
    if size(frameSet,2) ~= size(gazeTargets,2)
        error('estimateSceneParams:inputTypeError','The frameSet and gazeTargets must have equal second dimensions')
    end
    % Package the inputs as cell arrays
    videoStemName = {videoStemName};
    frameSet = {frameSet};
    gazeTargets = {gazeTargets};
    sceneArgs = {p.Results.sceneArgs};
end

% Define some constants
strategy = p.Results.searchStrategy;
verbose = p.Results.verbose;
nScenes = length(videoStemName);

% Set up empty cell arrays for the passed data if needed, and check that
% the lengths are correct. If these are empty, then the data will be loaded
% from files based upon the videoStemName
if isempty(p.Results.glintData)
    glintData = repmat({''},nScenes,1);
else
    glintData = p.Results.glintData;
end
if isempty(p.Results.perimeter)
    perimeter = repmat({''},nScenes,1);
else
    perimeter = p.Results.perimeter;
end
if isempty(p.Results.relativeCameraPosition)
    relativeCameraPosition = repmat({''},nScenes,1);
else
    relativeCameraPosition = p.Results.relativeCameraPosition;
end
if range([length(glintData),length(perimeter),length(relativeCameraPosition),nScenes])~=0
    error('estimateSceneParams:inputTypeError','One or more of the passed data cells (glint, perim, camera position) is not the right length')
end


%% Define model params
% This function has a dictionary of search parameters and stages. The
% key-value 'model' may be used to supply values that replace the defaults.
% This is typically done for x0 and bounds.
model = defineModelParams(nScenes, p.Results.model, p.Results.cameraTorsion, p.Results.cameraDepth, p.Results.corneaTorsion);


%% Create the scene objective functions
for ss = 1:nScenes
    
    % The arguments for createSceneGeometry for this video entry are a
    % combination of the eyeArgs, and the sceneArgs for this video
    if isempty(p.Results.sceneArgs)
        setupArgs = p.Results.eyeArgs;
    else
        setupArgs = [p.Results.eyeArgs,sceneArgs{ss}];
    end
    
    % Create the objective for this scene
    sceneObjects{ss} = sceneObj(...
        model, videoStemName{ss}, frameSet{ss}, gazeTargets{ss}, ...
        setupArgs, p.Results, 'verbose', verbose, ...
        'glintData',glintData{ss},'perimeter',perimeter{ss}, ...
        'relativeCameraPosition',relativeCameraPosition{ss});
    
end

% Set the initial value of x to x0
x = model.x0;


%% Define BADS search options
options = bads('defaults');          % Get a default OPTIONS struct
options.Display = 'off';             % Silence display output
options.UncertaintyHandling = 0;     % The objective is deterministic
options.MeshOverflowsWarning = Inf;  % Silences this warning
options.TolMesh = model.strategy.(strategy).TolMesh;


%% Set up the parallel pool
if p.Results.useParallel
    startParpool( p.Results.nWorkers, verbose );
end


%% Announce we are starting
ticObject = tic();
if verbose
    fprintf(['Estimating scene parameters. Started ' char(datetime('now')) '\n']);
end


%% Search across stages

% The number of stages
nStages = length(model.strategy.(strategy).stages);

% Looping over stages
for ii = 1:nStages
    
    % Announce
    if verbose
        str = sprintf('Starting stage %d of %d \n',ii,nStages);
        fprintf(str);
    end
    
    % Objective
    myObjAll = @(x) multiSceneObjective(x,sceneObjects,model,strategy,ii,p.Results.errorArgs,verbose);
    
    % Bounds
    [x,lb,ub,lbp,ubp] = setBounds(x,model,ii,strategy);
    
    % Search
    if any(lb~=ub)
        x = bads(myObjAll,x,lb,ub,lbp,ubp,model.func.nonbcon,options);
    end
    
    % Update the model objects with the final parameters
    myObjAll(x);
    
    % Plots
    if p.Results.savePlots
        fileNameSuffix = sprintf([p.Results.outputFileSuffix '_stage%02d'],ii);
        for ss = 1:nScenes
            sceneObjects{ss}.saveEyeModelMontage(fileNameSuffix,true,true,true);
            sceneObjects{ss}.saveModelFitPlot(fileNameSuffix);
            sceneObjects{ss}.saveRelCameraPosPlot(fileNameSuffix);
        end
    end
    
    % If instructed, use the fixation results to update the primary
    % position after each stage except the last
    if model.strategy.(strategy).useFixForPrimaryPos && ii<nStages
        poses = [];
        for ss = 1:nScenes
            fixationEyePose = sceneObjects{ss}.poseRegParams.t;
            poses = [poses fixationEyePose(1:2)'];
        end
        idx = model.func.fieldSetIdx('scene','primaryPosition');
        idx = model.scene.idxMultiScene(idx);
        x(idx) = poses;
    end
    
    
end


%% Ensure that the scene objects are defined using the final value of x
myObjAll(x);


%% Save the sceneGeometry and relative camera position
if p.Results.saveFiles
    for ss = 1:nScenes
        sceneObjects{ss}.saveSceneGeometry(p.Results.outputFileSuffix);
        sceneObjects{ss}.saveRelCameraPos(p.Results.outputFileSuffix);
    end
end

%% Alert the user that we are done with the routine
executionTime = toc(ticObject);
if verbose
    fprintf([num2str(executionTime) '\n']);
end


end


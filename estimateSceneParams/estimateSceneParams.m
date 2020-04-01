function estimateSceneParams(videoStemName, frameSet, gazeTargets, varargin)
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
%   The accuracy of the solution is improved by setting the cameraDepth
%   key-value. This specifies the distance of the camera from the eye (in
%   mm), as it is otherwise difficult for the solution to distinguish
%   between the depth position of the camera and properties of eye biometry
%   (rotation depth and corneal curvature). The routine includes a
%   regularization that penalizes the solution for departing from the x0
%   value for camera depth. The penalty rises as the square of the distance
%   in mm of the camera depth from the inital value, with the weight of the
%   penalty adjusted setting 'depthChangePenaltyWeight'. A value of zero
%   will remove the penalty.
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
%                           and optionally {'_headMotion.mat'}. The
%                           sceneGeometry output files will be saved at the
%                           same path location with the suffix:
%                               '_sceneGeometry.mat'
%   frameSet              - A 1xm vector, or a cell array of n such vectors
%                           (although m may differ in each vector). Each
%                           vector specifies the m frame indices (indexed
%                           from 1) which identify the set of frames from
%                           the scene to guide the search.
%   gazeTargets           - A 2xm matrix, or cell array of n such matrices,
%                           (although m may differ in each matrix). Each
%                           matrix provides the positions, in degrees of
%                           visual angle, of fixation targets that
%                           correspond to each of the frames. The visual
%                           angle of the stimuli should be adjusted for
%                           min/magnification produced by spectacle lenses
%                           worn by the subject prior to being passed to
%                           this routine.
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
%  'searchStrategy'       - Char vector. See the Description section for
%                           details, and the function defineModelParams.m
%  'cameraDepth'          - Scalar. Estimate of the distance of the nodal
%                           point of the camera from the corneal apex of
%                           the eye. Used to inform the x0 values for the
%                           search.
%  'depthChangePenaltyWeight' - Scalar. Adjustment of camera depth away
%                           from the x0 values is penalized. The effect of
%                           the weight (w) is such that an x% change in
%                           camera depth produces a [(w*x)^2]% increase in
%                           the objective. A weight value of zero removes
%                           the regularization.
%  'model'                  Structure. Setting fields of this structure to
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
%                               { {'spectacleLens',-3} ,{'spectacleLens',-3} }
%  'errorArgs'            - Cell array. These are key-value pairs that are
%                           passed to method 'updateError' of sceneObj.
%                           Typically, one might pass pre-computed values
%                           of poseRegParams or vecRegParams, or adjust the
%                           weighting function for the types of errors.
%  'useFixForPrimaryPos'  - Logical. The primary position of the eye is
%                           specified as part of the eye model, and
%                           controls the pseudo-torsion of the eye
%                           following Listing's Law. Because eye rotations
%                           are defined with respect to the alignment of
%                           the optical axis of the camera with the optical
%                           axis of the eye, the primary position of the
%                           eye in the scene coordinate frame will vary
%                           with camera position. The routine is able to
%                           search across primary position values for each
%                           scene. This flag causes the model to be updated
%                           after each stage such that the primary position
%                           for the eye in a given scene is set equal to
%                           the modeled position of the eye when it is
%                           fixating position [0 0] of a fixation array.
%                           The use of this flag is justified when the
%                           subject is allowed to adjust their head to most
%                           comfortably view the center of a fixation array
%                           prior to recording, presumably placing the
%                           center of fixation at their primary position.
%  'multiSceneNorm'       - Scalar. Defines the metric used to combine
%                           errors across scenes. Defaults to a value of 1
%                           and thus the L1 norm.
%  'TolMesh'              - Scalar. The precision with which the parameters
%                           are adjusted within the BADS search. A value of
%                           1e-2 provides a good trade-off between search
%                           time and useful precision in the eye and scene
%                           parameters.
%
% Outputs
%	None. Although sceneGeometry files and diagnostic plots are saved for
%	each scene at the location specified by the videoStemName paths.
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
    sourceSceneName = 'GazeCal01_sceneGeometry.mat';
    load(sourceSceneName,'sceneGeometry');
    model.eye.x0 = sceneGeometry.meta.estimateSceneParams.x(1:7);
    eyeArgs = sceneGeometry.meta.estimateSceneParams.p.eyeArgs;
    errorArgs = {...
        'errorReg',[1 1 0 0], ...
        'poseRegParams',sceneGeometry.meta.estimateSceneParams.poseRegParams,...
        'vectorRegParams',sceneGeometry.meta.estimateSceneParams.vectorRegParams};

    % This is the video for which we wish to create a sceneGeometry file. 
    % We select frames to guide the search, using a fixation period before
    % the scan and a distributed set of gaze positions after the scan start
    videoStemName = 'tfMRI_MOVIE_PA_run03';
    [frameSet1, gazeTargets1] = selectFrames.gazePre(videoStemName);
    [frameSet2, gazeTargets2] = selectFrames.grid(videoStemName);
    frameSet = [frameSet1 frameSet2];
    gazeTargets = [gazeTargets1 gazeTargets2];

    % Perform the search
    estimateSceneParams(videoStemName, frameSet, gazeTargets, ...
        'searchStrategy','sceneSync','cameraDepth',140,'model',model,...
        'eyeArgs',eyeArgs,'sceneArgs',sceneArgs,'errorArgs',errorArgs);
%}



%% input parser
p = inputParser; p.KeepUnmatched = true;

% Required
p.addRequired('videoStemName',@(x)(iscell(x) || ischar(x)));
p.addRequired('frameSet',@(x)(iscell(x) || isnumeric(x)));
p.addRequired('gazeTargets',@(x)(iscell(x) || isnumeric(x)));

% Optional display and I/O params
p.addParameter('verbose',true,@islogical);

% Optional flow control params
p.addParameter('useParallel',true,@islogical);
p.addParameter('nWorkers',[],@(x)(isempty(x) || isnumeric(x)));

% Optional environment parameters
p.addParameter('tbSnapshot',[],@(x)(isempty(x) || isstruct(x)));
p.addParameter('timestamp',char(datetime('now')),@ischar);
p.addParameter('username',char(java.lang.System.getProperty('user.name')),@ischar);
p.addParameter('hostname',char(java.net.InetAddress.getLocalHost.getHostName),@ischar);

% Optional analysis params
p.addParameter('searchStrategy','gazeCal',@ischar);
p.addParameter('cameraDepth',120,@isnumeric);
p.addParameter('depthChangePenaltyWeight',0.5,@isscalar);
p.addParameter('model',[],@isstruct);
p.addParameter('eyeArgs',{},@iscell);
p.addParameter('sceneArgs',{},@iscell);
p.addParameter('errorArgs',{},@iscell);
p.addParameter('useFixForPrimaryPos',true,@islogical);
p.addParameter('multiSceneNorm',1,@isscalar);
p.addParameter('TolMesh',1e-2,@isscalar);

% parse
p.parse(videoStemName, frameSet, gazeTargets, varargin{:})


%% Check inputs

% Catch some possible input errors
if iscell(videoStemName)
    % If the videoNameStem is a cell array, make sure that all inputs are cell
    % arrays
    if ~iscell(frameSet) || ~iscell(frameSet)
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
end

% Pull these out for code legibility
verbose = p.Results.verbose;
nScenes = length(videoStemName);


%% Define model params
% This function has the dictionary of search parameters and stages.
% The key-value 'model' may be used to supply values that replace the
% defaults. This is typically done for x0 and bounds.
model = defineModelParams(nScenes, p.Results.model, p.Results.cameraDepth);


% Loop through the scenes and create scene objective functions
for ss = 1:nScenes
    
    % The arguments for createSceneGeometry for this video entry are a
    % combination of the eyeArgs, and the sceneArgs for this video
    if isempty(p.Results.sceneArgs)
        setupArgs = p.Results.eyeArgs;
    else
        setupArgs = [p.Results.eyeArgs,p.Results.sceneArgs{ss}];
    end
    
    % Create the objective for this scene
    mySceneObjects{ss} = sceneObj(...
        model, videoStemName{ss}, frameSet{ss}, gazeTargets{ss}, ...
        setupArgs, p.Results.errorArgs, p.Results, ...
        'verbose', verbose);
    
end

% Set the initial value of x to x0
x = model.x0;


%% Anonymous functions for the search
% An objective function which is the norm of all objective functions

myObjAll = @(x) multiSceneObjective(x,mySceneObjects,model,p.Results.multiSceneNorm,verbose);


%% Define BADS search options
options = bads('defaults');          % Get a default OPTIONS struct
options.Display = 'off';             % Silence display output
options.UncertaintyHandling = 0;     % The objective is deterministic
options.TolMesh = p.Results.TolMesh; % Typically 1e-3 is plenty precise


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

% The number of sta
nStages = length(model.stages.(p.Results.searchStrategy));
for ii = 1:nStages
    
    % Announce
    if verbose
        str = sprintf(['Starting stage %d of %d \n'],ii,nStages);
        fprintf(str);
    end
    
    % Bounds
    [x,lb,ub,lbp,ubp] = setBounds(x,model,ii,p.Results.searchStrategy);
    
    % Search
    x = bads(myObjAll,x,lb,ub,lbp,ubp,model.func.nonbcon,options);
    
    % Plots
    fileNameSuffix = sprintf('_stage%02d',ii);
    for ss = 1:nScenes
        mySceneObjects{ss}.saveEyeModelMontage(fileNameSuffix);
        mySceneObjects{ss}.saveModelFitPlot(fileNameSuffix);
    end        
    
    % If instructed, use the fixation results to update the primary
    % position after each stage except the last
    if p.Results.useFixForPrimaryPos && ii<nStages
        poses = [];
        for ss = 1:nScenes
            fixationEyePose = mySceneObjects{ss}.fixationEyePose;
            poses = [poses fixationEyePose(1:2)'];
        end
        idx = model.func.fieldSetIdx('scene','primaryPosition');
        idx = model.scene.idxMultiScene(idx);
        x(idx) = poses;
    end

    
end


%% Ensure that the scene objects are defined using the final value of x
myObjAll(x);


%% Save the sceneGeometry
for ss = 1:nScenes
    mySceneObjects{ss}.saveSceneGeometry('');
end


%% Alert the user that we are done with the routine
executionTime = toc(ticObject);
if verbose
    fprintf([num2str(executionTime) '\n']);
end


end


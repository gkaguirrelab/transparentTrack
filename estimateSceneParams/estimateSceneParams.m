function estimateSceneParams(videoStemName, frameSet, gazeTargets, varargin)
% Estimate eye and scene geometry parameters for one or more scenes
%
% Syntax:
%  sceneGeometry = estimateSceneParams(pupilFileName, perimeterFileName, glintFileName, sceneGeometryFileName, varargin)
%
% Description:
%   The appearance of the eye in a camera image is influenced by several
%   parameters, including biometric properties of the eye and the position
%   of the camera. This function uses a set of observations of the pupil
%   and the glint to estimate these scene parameters, including the location of visual targets
%   that the eye was fixated upon during each of the observations. A set of
%   measurements made with the camera in a fixed location constitutes a
%   "scene". The routine may be supplied with information from more than
%   one scene. The routine forces the eye biometric parameters to be the
%   same in all scenes, but allows the scene parameters (camera position)
%   to vary between scenes.
%
%   The accuracy of the solution is greatly improved by supplying an
%   accurate initial value for the distance of the camera from the eye in
%   the sceneParamsX0 vector. Unless the scenes feature the camera in very
%   different positions, it is difficult for the solution to distinguish
%   between the depth position of the camera and properties of eye biometry
%   (rotation depth and corneal curvature). The routine includes a
%   regularization that penalizes the solution for departing from the x0
%   value for camera depth. The penalty rises as the square of the distance
%   in mm of the camera depth from the inital value, with the weight of the
%   penalty adjusted using the key 'depthChangePenaltyWeight'. Setting this
%   weight to zero will remove the penalty.
%
% Inputs:
%	videoStemName         - Cell array of length n of char vectors. Full
%                           path to the n video files from which the scene
%                           observations have been derived. The stem name
%                           should omit the "_gray.avi" suffix that is
%                           usually present in the names of these video
%                           files. The routine assumes that files to be
%                           used in the analysis have the suffixes:
%                               {'_correctedPerimeter.mat','_glint.mat', ...
%                               '_pupil.mat'}
%                           The sceneGeometry output files will be saved at
%                           the same path location with the suffix:
%                               '_sceneGeometry.mat'
%   frameSet              - Cell array of length n of numeric vectors. Each
%                           vector specifies the m frame indices (indexed
%                           from 1) which identify the set of frames from
%                           the scene to guide the search.
%   gazeTargets           - Cell array of length n of  2xm matrices that
%                           provide the positions, in degrees of visual
%                           angle, of fixation targets that correspond to
%                           each of the frames. The visual angle of the
%                           stimuli should be adjusted for
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
%  'eyeArgs'              - Cell array. This is a cell array of key-value
%                           pairs to be used when generating the eye model
%                           for this subject. An example input is:
%                           {'axialLength',23.45,'sphericalAmetropia',-0.5}
%                           Note that these key values:
%                               {'kvals','rotationCenterScalers',...
%                               'primaryPosition'}
%                           should not be passed, as they will be
%                           supplanted parameters in the search. If you
%                           wish to provide and lock these parameters in
%                           the search, do so by providing the values in
%                           the x0 arrays, and setting the corresponding
%                           values in the param bounds to zero.
%  'sceneArgs'            - Cell array of length n of cell arrays. Each
%                           array contians key-value pairs to be used when
%                           generating the scene model for a particular
%                           video. An example input (for 2 scenes ) is:
%                           { {'spectacleLens',-3} , {'spectacleLens',-3} }
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
%  'eyeParamsX0'          - 1x7 vector. Provides the initial values for the
%                           parameters that control the biometric
%                           properties of the eye The order of these
%                           parameters is specified in the hard-coded
%                           variable: eyeParamLabels.
%  'eyeParamsBounds'      - 1x7 vector. These values are added to and
%                           subtracted from the x0 vector to define the
%                           hard limits of the search.
%  'sceneParamsX0'        - Cell array of length n of 1x6 vectors. Provides
%                           the initial values for the paramters that
%                           control camera position in each scene. It is
%                           reasonable to set all of these parameter values
%                           to zero, with the exception of the depth
%                           parameter, which should be as accurate as
%                           possible. The order of the parameters is
%                           specified in the variable: sceneParamLabels
%  'sceneParamsBounds'    - 1x6 vector. These values are added to and
%                           subtracted from the x0 vector to define the
%                           hard limits of the search. The values in this
%                           vector are used for all scenes.
%  'eyePoseLB/UB'         - 1x4 vector. Upper / lower bounds on the eyePose
%                           [azimuth, elevation, torsion, pupil radius].
%                           Biological limits in eye rotation and pupil
%                           size would suggest boundaries of [±35, ±25, 0,
%                           0.25-4]. Note, however, that these angles are
%                           relative to the center of projection, not the
%                           primary position of the eye. Therefore, in
%                           circumstances in which the camera is viewing
%                           the eye from an off-center angle, the bounds
%                           will need to be shifted accordingly.
%  'errorReg'             - 1x4 vector. Passed to calcGlintGazeError and
%                           determines the relative weighting of the
%                           different types of error in scene prediction.
%  'multiSceneNorm'       - Scalar. Defines the metric used to combine
%                           errors across scenes. Defaults to a value of 1
%                           and thus the L1 norm.
%  'TolMesh'              - Scalar. The precision with which the parameters
%                           are adjusted within the BADS search. A value of
%                           1e-2 provides a good trade-off between search
%                           time and useful precision in the eye and scene
%                           parameters.
%  'depthChangePenaltyWeight' - Scalar. Adjustment of camera depth away
%                           from the x0 values is penalized. The effect of
%                           the weight (w) is such that an x% change in
%                           camera depth produces a [(w*x)^2]% increase in
%                           the objective. A weight value of zero removes
%                           the regularization.
%
% Outputs
%	None. Although sceneGeometry files and diagnostic plots are saved for
%	each scene at the location specified by the videoStemName paths.
%
% Examples:
%{
    % Get the DropBox base directory
    dropboxBaseDir = getpref('eyeTrackTOMEAnalysis','dropboxBaseDir');
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
    sceneArgs = {'','','','',''};
    sceneParamsX0 = {...
        [0         0  -11.6168   -2.5423    2.4580  140], ...
        [0         0   -4.7981   -2.2882    1.4291  140], ...
        [0         0    0.3887   -2.1856    1.0091  140], ...
        [0         0   10.2627   -3.3679    3.5300  140], ...
        [0         0    7.7354   -0.8179   -2.7592  140]};
    eyeArgs = {'axialLength',23.45,'sphericalAmetropia',-0.5};
    eyeParamsX0 = [41.8000   42.8000         0         2.5         0    0.9357    0.9575];
    estimateSceneParams(videoStemName, frameSet, gazeTargets, ...
        'sceneArgs',sceneArgs, ...
        'sceneParamsX0', sceneParamsX0, ...
        'eyeArgs',eyeArgs, ...
        'eyeParamsX0',eyeParamsX0);
%}



%% input parser
p = inputParser; p.KeepUnmatched = true;

% Required
p.addRequired('videoStemName',@iscell);
p.addRequired('frameSet',@iscell);
p.addRequired('gazeTargets',@iscell);

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
p.addParameter('searchStrategy','twoStage',@ischar);
p.addParameter('eyeArgs',{''},@iscell);
p.addParameter('sceneArgs',{''},@iscell);
p.addParameter('keyArgs',{''},@iscell);
p.addParameter('useFixForPrimaryPos',true,@islogical);
p.addParameter('relCamPosParamsX0',[0, 0, 0, 0],@isnumeric);
p.addParameter('relCamPosParamsBounds',[30, 30, 30, 30],@isnumeric);
p.addParameter('eyeParamsX0',[44.2410, 45.6302, 0, 2.5000, 0, 1, 1],@isnumeric);
p.addParameter('eyeParamsBounds',[5, 5, 90, 5, 5, 0.25, 0.15],@isnumeric);
p.addParameter('sceneParamsX0',{[0 0 0 0 0 120]},@iscell);
p.addParameter('sceneParamsBounds',[10 10 20 20 20 20],@isnumeric);
p.addParameter('eyePoseLB',[-89,-89,0,0.1],@isnumeric);
p.addParameter('eyePoseUB',[89,89,0,4],@isnumeric);
p.addParameter('errorReg',[1 2 4 2],@isnumeric);
p.addParameter('multiSceneNorm',1,@isscalar);
p.addParameter('TolMesh',1e-2,@isscalar);
p.addParameter('depthChangePenaltyWeight',0.5,@isscalar);

% parse
p.parse(videoStemName, frameSet, gazeTargets, varargin{:})



%% Define model parameters

% The number of scenes to be fit
nScenes = length(videoStemName);

% The relative camera position parameters to be adjusted
relCamPosParamLabels = {'scanner_timeShift','scanner_azi','scanner_ele','scanner_torsion'};
nRelCamPosParams = length(relCamPosParamLabels);

% The eye model parameters that are to be adjusted
eyeParamLabels = {'cornea_K1','cornea_K2','cornea_torsion','cornea_tilt','cornea_tip','rotationDepth_joint','rotationDepth_diff'};
nEyeParams = length(eyeParamLabels);

% The scene parameters that are to be adjusted
sceneParamLabels = {'primaryPosition_azi','primaryPosition_ele','camera_torsion','camera_horizontal','camera_vertical','camera_depth'};
nSceneParams = length(sceneParamLabels);

% The eye parameters are shared by all scenes, and each scene gets its own
% set of scene parameters, yielding this many total parameters in the
% optimization.
nTotalParams = nRelCamPosParams + nEyeParams + nSceneParams * nScenes;


%% Define parameter search sets
% These are vectors that are used to lock or free different sets of
% parameters in the search

% The searchSet vectors will be of the same size as the parameter vector
blankSearch = uint8(zeros(1,nTotalParams));

% An anonymous function that expands the input index vector [a, b, ...]
% into a vector given k eye+relCamPosParams, s sceneParams, and n scenes:
%	[ e+a+(s*0), e+b+(s*0), e+a+(s*1), e+b+(s*1), ... e+a+(s*(n-1)), e+b+(s*(n-1)) ]
sceneIdxRep = @(idx) nRelCamPosParams + nEyeParams + repmat((0:nScenes-1)*nSceneParams,1,length(idx)) + ...
    cell2mat(arrayfun(@(x) repmat(x,1,nScenes),idx,'UniformOutput',false));

% The search vectors
scannerTimeSet = blankSearch; scannerTimeSet(1)=1;
scannerRotationSet = blankSearch; scannerRotationSet(2:4)=1;
corneaSet = blankSearch; corneaSet(5:9)=1;
rotationSet = blankSearch; rotationSet(10:11)=1;
primaryPosSet = blankSearch; primaryPosSet(sceneIdxRep(1:2)) = 1;
cameraTorsionSet = blankSearch; cameraTorsionSet(sceneIdxRep(3)) = 1;
cameraPlaneTransSet = blankSearch; cameraPlaneTransSet(sceneIdxRep(4:5)) = 1;
cameraDepthTransSet = blankSearch; cameraDepthTransSet(sceneIdxRep(6)) = 1;

% Create stages of search sets, based upon the identified strategy
switch p.Results.searchStrategy
    case 'twoStage'
        searchSets = {...
            logical(rotationSet + cameraTorsionSet + cameraPlaneTransSet + cameraDepthTransSet), ...
            logical(corneaSet + rotationSet + primaryPosSet + cameraTorsionSet + cameraPlaneTransSet + cameraDepthTransSet) ...
            };
        searchSetLabels = {'Camera position and eye rotation','All scene and eye parameters'};

    case 'sceneSync'
        searchSets = {...
            logical( cameraTorsionSet + cameraPlaneTransSet + cameraDepthTransSet), ...
            logical( scannerRotationSet), ...
            logical( scannerTimeSet + primaryPosSet) ...
            };
        searchSetLabels = {'Camera position','Scanner position','Scanne time and Primary eye position'};
    otherwise
        error('Not a recognized searchStrategy');
end


nStages = length(searchSets);


%% Set up the search options and objectives

% The eye parameters are the initial entries in x0 and bounds
x0 = [p.Results.relCamPosParamsX0, p.Results.eyeParamsX0];
xBounds = [p.Results.relCamPosParamsBounds, p.Results.eyeParamsBounds];

% These key values are passed to calcGlintGazeError
keyVals = {...
    'eyePoseLB', p.Results.eyePoseLB,...
    'eyePoseUB', p.Results.eyePoseUB,...
    'errorReg',p.Results.errorReg, ...
    p.Results.keyArgs{:} ...
    };

% Loop through the scenes and create scene objective functions
for ss = 1:nScenes
    
    % The arguments for createSceneGeometry for this video entry are a
    % combination of the eyeArgs, and the sceneArgs for this video
    setupArgs = [p.Results.eyeArgs,p.Results.sceneArgs{ss}];
    
    % Create the objective for this scene
    mySceneObjects{ss} = sceneObj(...
        videoStemName{ss}, frameSet{ss}, gazeTargets{ss}, ...
        setupArgs, keyVals, p.Results, ...
        'verbose', p.Results.verbose);
    
    % Assemble the x0 by concatenating the x0 scene params
    x0 = [x0, p.Results.sceneParamsX0{ss}];
    
    % Assemble xBounds by concatenating the scene bounds
    xBounds = [xBounds, p.Results.sceneParamsBounds];
end

% Set the initial value of x to x0
x = x0;


%% Anonymous functions for the search

% A regularization function that penalizes changes in depth from the x0
% values
depthChangePenaltyMultiplier = @(x) (1 + p.Results.depthChangePenaltyWeight * norm( (x(logical(cameraDepthTransSet)) - x0(logical(cameraDepthTransSet))) ./ x0(logical(cameraDepthTransSet)) ))^2;

% An objective function which is the norm of all objective functions
myObjAll = @(x) multiSceneObjective(x,mySceneObjects,nEyeParams,nRelCamPosParams,nSceneParams,depthChangePenaltyMultiplier,p.Results.multiSceneNorm,p.Results.verbose);

% A non-linear constraint for the BADS search that requires first value of
% the corneal curvature (K1) to be less than the second value (K2) Note
% that NONBCON takes a matrix input, which is why we perform this
% calculation over the first dimension.
nonbcon = @(x) x(:,1) > x(:,2);


%% Define BADS search options
options = bads('defaults');          % Get a default OPTIONS struct
options.Display = 'off';             % Silence display output
options.UncertaintyHandling = 0;     % The objective is deterministic
options.TolMesh = p.Results.TolMesh; % Typically 1e-3 is plenty precise


%% Set up the parallel pool
if p.Results.useParallel
    startParpool( p.Results.nWorkers, p.Results.verbose );
end


%% Announce we are starting
ticObject = tic();
if p.Results.verbose
    fprintf(['Estimating scene parameters. Started ' char(datetime('now')) '\n']);
end


%% Search across stages
for ii = 1:nStages
    
    % Announce
    if p.Results.verbose
        str = sprintf(['Starting stage %d of %d, ' searchSetLabels{ii} '\n'],ii,nStages);
        fprintf(str);
    end
    
    % Bounds
    [x,lb,ub,lbp,ubp] = setBounds(x,xBounds,searchSets{ii});
    
    % Search
    x = bads(myObjAll,x,lb,ub,lbp,ubp,nonbcon,options);
    
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
        x(find(primaryPosSet)) = poses;
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
if p.Results.verbose
    fprintf([num2str(executionTime) '\n']);
end


end




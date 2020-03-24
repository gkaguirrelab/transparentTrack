function estimateSceneParams(videoStemName, frameSet, gazeTargets, varargin)
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
%   that the eye was fixated upon during each of the observations. It is
%   also valuable to provide a reasonably accurate value for the distance
%   of the camera from the eye in the sceneParamsX0 vector.
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
%  'saveDiagnosticPlot'   - Logical.
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
%                           the glint and pupil perimeter. If needed, the
%                           visual angle of the stimuli will be adjusted
%                           for min/magnification produced by spectacle
%                           lenses worn by the subject.
%  'fixSpectacleLens'     - Scalar. This parameter is used to handle the
%                           unusual circumstance in which the eye viewing a
%                           fixation array is behind a spectacle lens, but
%                           the eye being modeled for pupil appearance is
%                           not. Setting this parameter causes the routine
%                           to calculate a magnification factor for the
%                           fixation target array, but does not apply this
%                           spectacle to the eye being modeled.
%  'searchThresh'         - Scalar. Search iterations will be repeated
%                           until the fVal is below this threshold, or the
%                           maximum number of search iterations has been
%                           reached.
%  'searchIterations'     - Scalar integer. The maximum number of search
%                           iterations to conduct.
%  'sceneParamsX0'        - 1x11 vector. Scene parameters to use as the
%                           starting point for the search.
%  'sceneParamsToLock'    - 1x11 vector. Any param positions set to unity
%                           will be locked at their x0 values.
%
% Outputs
%	sceneGeometry         - A structure that contains the components of the
%                           projection model.
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
        [0.00, 0.00, -6.88, -2.59, 2.21, 140.81], ...
        [0.00, 0.00, 2.36, -2.35, 1.19, 140.83], ...
        [0.00, 0.00, 0.00, 0.00, 0.00, 140.00], ...
        [0.00, 0.00, 0.00, 0.00, 0.00, 140.00], ...
        [0.00, 0.00, 0.00, 0.00, 0.00, 140.00]};
    eyeArgs = {'axialLength',23.45,'sphericalAmetropia',-0.5};
    eyeParamsX0 = [41.80, 42.80, 0.00, 0.00, 0.00, 0.91, 0.97];
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
p.addParameter('saveDiagnosticPlot',true,@islogical);

% Optional flow control params
p.addParameter('useParallel',true,@islogical);
p.addParameter('nWorkers',[],@(x)(isempty(x) || isnumeric(x)));

% Optional environment parameters
p.addParameter('tbSnapshot',[],@(x)(isempty(x) || isstruct(x)));
p.addParameter('timestamp',char(datetime('now')),@ischar);
p.addParameter('username',char(java.lang.System.getProperty('user.name')),@ischar);
p.addParameter('hostname',char(java.net.InetAddress.getLocalHost.getHostName),@ischar);

% Optional analysis params
p.addParameter('eyeArgs',{''},@iscell);
p.addParameter('sceneArgs',{''},@iscell);
p.addParameter('searchThresh',1.0,@isscalar);
p.addParameter('useFixForPrimaryPos',true,@islogical);
p.addParameter('eyeParamsX0',[44.2410, 45.6302, 0, 2.5000, 0, 1, 1],@isnumeric);
p.addParameter('eyeParamsBounds',[5, 5, 90, 5, 5, 0.25, 0.15],@isnumeric);
p.addParameter('sceneParamsX0',{[0 0 0 0 0 120]},@iscell);
p.addParameter('sceneParamsBounds',[10 10 20 20 20 20],@isnumeric);
p.addParameter('eyePoseLB',[-89,-89,0,0.1],@isnumeric);
p.addParameter('eyePoseUB',[89,89,0,4],@isnumeric);
p.addParameter('errorReg',[1 2 4 2],@isnumeric);
p.addParameter('multiSceneNorm',1,@isscalar);


% parse
p.parse(videoStemName, frameSet, gazeTargets, varargin{:})



%% Define model parameters

nScenes = length(videoStemName);

eyeParamLabels = {'cornea_K1','cornea_K2','cornea_torsion','cornea_tilt','cornea_tip','rotationDepth_joint','rotationDepth_diff'};
nEyeParams = length(eyeParamLabels);

sceneParamLabels = {'primaryPosition_azi','primaryPosition_ele','camera_torsion','camera_horizontal','camera_vertical','camera_depth'};
nSceneParams = length(sceneParamLabels);

nTotalParams = nEyeParams + nSceneParams * nScenes;

%% Define parameter search sets
blankSearch = zeros(1,nTotalParams);
sceneIdxRep = @(idx) nEyeParams + repmat((0:nScenes-1)*nSceneParams,1,length(idx)) + ...
    cell2mat(arrayfun(@(x) repmat(x,1,nScenes),idx,'UniformOutput',false));

corneaSet = blankSearch; corneaSet(1:5)=1;
rotationSet = blankSearch; rotationSet(6:7)=1;
primaryPosSet = blankSearch; primaryPosSet(sceneIdxRep(1:2)) = 1;
cameraTorsionSet = blankSearch; cameraTorsionSet(sceneIdxRep(3)) = 1;
cameraPlaneTransSet = blankSearch; cameraPlaneTransSet(sceneIdxRep(4:5)) = 1;
cameraDepthTransSet = blankSearch; cameraDepthTransSet(sceneIdxRep(6)) = 1;


%% Set up the search components

% The eye parameters are the initial entries in x0 and bounds
x0 = p.Results.eyeParamsX0;
xBounds = p.Results.eyeParamsBounds;

% These key values are passed to calcGlintGazeError
keyVals = {...
    'eyePoseLB', p.Results.eyePoseLB,...
    'eyePoseUB', p.Results.eyePoseUB,...
    'errorReg',p.Results.errorReg ...
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

% Create the objective function which is the norm of all objective
% functions
myObjAll = @(x) multiSceneObjective(x,mySceneObjects,nEyeParams,nSceneParams,p.Results.multiSceneNorm,p.Results.verbose);


%% Define BADS search options
options = bads('defaults');          % Get a default OPTIONS struct
options.Display = 'off';             % Silence display output
options.UncertaintyHandling = 0;     % The objective is deterministic

% Define a non-linear constraint for the BADS search that requires first
% value of the corneal curvature (K1) to be less than the second value (K2)
% Note that NONBCON expects a matrix input.
nonbcon = @(x) x(:,1) > x(:,2);


%% Announce we are starting
ticObject = tic();
if p.Results.verbose
    fprintf(['Estimating scene parameters. Started ' char(datetime('now')) '\n']);
end


%% Set up the parallel pool
if p.Results.useParallel
    startParpool( p.Results.nWorkers, p.Results.verbose );
end


%% Search
x = x0;

% Define the initial search set
searchSet = rotationSet + cameraTorsionSet + cameraPlaneTransSet + cameraDepthTransSet;

% Set the bounds
[x,lb,ub,lbp,ubp] = setBounds(x,xBounds,searchSet);

% Perform the initial search
x = bads(myObjAll,x,lb,ub,lbp,ubp,nonbcon,options);

% Save out diagnostic plots for this stage
fileNameSuffix = '_stage1';
for ss = 1:nScenes
    mySceneObjects{ss}.saveEyeModelMontage(fileNameSuffix);
    mySceneObjects{ss}.saveModelFitPlot(fileNameSuffix);
end

% Use the fixation results to define eye primary position
if p.Results.useFixForPrimaryPos
    poses = [];
    for ss = 1:nScenes
        fixationEyePose = mySceneObjects{ss}.fixationEyePose;
        poses = [poses fixationEyePose(1:2)'];
    end
    x(find(primaryPosSet)) = poses;
end

% Define the full search set
searchSet = corneaSet + rotationSet + primaryPosSet + cameraTorsionSet + ... 
    cameraPlaneTransSet + cameraDepthTransSet;

% Set the bounds
[x,lb,ub,lbp,ubp] = setBounds(x,xBounds,searchSet);

% Perform the search
x = bads(myObjAll,x,lb,ub,lbp,ubp,nonbcon,options);
myObjAll(x);

% Save the sceneGeometry and plots
fileNameSuffix = '_stage2';
for ss = 1:nScenes
    mySceneObjects{ss}.saveEyeModelMontage(fileNameSuffix);
    mySceneObjects{ss}.saveModelFitPlot(fileNameSuffix);
    mySceneObjects{ss}.saveSceneGeometry('');
end

% Get the execution time
executionTime = toc(ticObject);

%% alert the user that we are done with the routine
if p.Results.verbose
    fprintf([num2str(executionTime) '\n']);
end

end






%%%%%%%%%%%% LOCAL FUNCTIONS



function [x,lb,ub,lbp,ubp] = setBounds(x,bb,searchSet)

% Apply the bounds
lb = x - bb;
lbp = x - bb./2;
ubp = x + bb./2;
ub = x + bb;

% Lock params
lb(~searchSet) = x(~searchSet);
ub(~searchSet) = x(~searchSet);
lbp(~searchSet) = x(~searchSet);
ubp(~searchSet) = x(~searchSet);

% Keep the corneal axis (torsion) between +-90
lb(3) = max([lb(3) -90]);
ub(3) = min([ub(3) 90]);
lbp(3) = max([lbp(3) -90]);
ubp(3) = min([ubp(3) 90]);

% Ensure that x is within bounds
x=min([ub; x]);
x=max([lb; x]);

end






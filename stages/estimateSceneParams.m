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
    sceneGeometryFileName = '/Users/aguirre/Dropbox (Aguirre-Brainard Lab)/TOME_processing/session2_spatialStimuli/TOME_3015/032417/EyeTracking/GazeCal03_sceneGeometry.mat';
    gazeTargets = [ -7, 0, -7, 7, 7, 0, -7, 0, 7 ; 0, -7, -7, 0, -7, 7, 7, 0, 7];
    frameSet = [ 679, 884, 1180, 1250, 1571, 1663, 1809, 2004, 2075 ];
    varargin = {'axialLength',23.45,'sphericalAmetropia',-0.5,'measuredCornealCurvature',[41.80,42.80,178]};
    estimateSceneParams(pupilFileName, perimeterFileName, glintFileName, sceneGeometryFileName, 'frameSet', frameSet, 'gazeTargets', gazeTargets, varargin{:});
%}
%{
    perimeterFileName = '/Users/aguirre/Dropbox (Aguirre-Brainard Lab)/TOME_processing/session2_spatialStimuli/TOME_3021/060917/EyeTracking/GazeCal02_correctedPerimeter.mat';
    glintFileName = '/Users/aguirre/Dropbox (Aguirre-Brainard Lab)/TOME_processing/session2_spatialStimuli/TOME_3021/060917/EyeTracking/GazeCal02_glint.mat';
    pupilFileName = '/Users/aguirre/Dropbox (Aguirre-Brainard Lab)/TOME_processing/session2_spatialStimuli/TOME_3021/060917/EyeTracking/GazeCal02_pupil.mat';
    sceneGeometryFileName = '/Users/aguirre/Dropbox (Aguirre-Brainard Lab)/TOME_processing/session2_spatialStimuli/TOME_3021/060917/EyeTracking/GazeCal02_sceneGeometry.mat';
    gazeTargets = [ -7, -7, 7, 0, 0, 7, 0, 7, -7 ; 0, 7, -7, -7, 0, 7, 7, 0, -7];
    frameSet = [ 730, 882, 971, 1114, 1250, 1382, 1467, 1593, 1672 ];
    varargin = {'axialLength',25.29,'sphericalAmetropia',-5.25,'contactLens',-5.25,'sceneParamsX0',[ 10.2446 -3.2242 -3.8030 124.1518 0.8628 0.9635 1.0488 0.9931 29.3357 ]};
    estimateSceneParams(pupilFileName, perimeterFileName, glintFileName, sceneGeometryFileName, 'frameSet', frameSet, 'gazeTargets', gazeTargets, varargin{:});
%}

%% input parser
p = inputParser; p.KeepUnmatched = true;

% Required
p.addRequired('pupilFileName',@ischar);
p.addRequired('perimeterFileName',@ischar);
p.addRequired('glintFileName',@ischar);
p.addRequired('sceneGeometryFileName',@ischar);

% Optional display and I/O params
p.addParameter('verbose',true,@islogical);
p.addParameter('grayVideoName','',@(x)(isempty(x) | ischar(x)));
p.addParameter('pupilFileToVideoSuffixSwitch',{'_pupil.mat','_gray.avi'},@iscell);

% Optional flow control params
p.addParameter('useParallel',true,@islogical);
p.addParameter('nWorkers',[],@(x)(isempty(x) || isnumeric(x)));

% Optional environment parameters
p.addParameter('tbSnapshot',[],@(x)(isempty(x) || isstruct(x)));
p.addParameter('timestamp',char(datetime('now')),@ischar);
p.addParameter('username',char(java.lang.System.getProperty('user.name')),@ischar);
p.addParameter('hostname',char(java.net.InetAddress.getLocalHost.getHostName),@ischar);

% Optional analysis params
p.addParameter('searchIterations',2,@isscalar);
p.addParameter('sceneParamsX0',[0 0 0 120 1 1 1 1 0],@isnumeric);
p.addParameter('lockDepth',false,@islogical);
p.addParameter('eyePoseLB',[-89,-89,0,0.1],@isnumeric);
p.addParameter('eyePoseUB',[89,89,0,4],@isnumeric);
p.addParameter('frameSet',[],@(x)(isempty(x) | isnumeric(x)));
p.addParameter('gazeTargets',[],@(x)(isempty(x) | isnumeric(x)));
p.addParameter('fixSpectacleLens',[],@(x)(isempty(x) | isnumeric(x)));

% parse
p.parse(pupilFileName, perimeterFileName, glintFileName, sceneGeometryFileName, varargin{:})



%% Error if gazeTargets or frameSet are empty
% Could add code here to derive these automatically if not set
if isempty(p.Results.frameSet) || isempty(p.Results.gazeTargets)
    warning('estimateSceneParams:undefinedFramesOrTargets','The routine currently requires that the frameSet and gazeTargets key-values be defined; returning');
    return
end


%% Announce we are starting
ticObject = tic();
if p.Results.verbose
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
% lens also has a (smaller) magnification effect, the lens rotates with the
% eye. Thus, eye movement amplitude is not altered.
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

frameSet = p.Results.frameSet;

% Extract the frames we want
perimeter.data = perimeter.data(frameSet);
ellipseRMSE = pupilData.initial.ellipses.RMSE(frameSet);
glintData.X = glintData.X(frameSet); glintData.Y = glintData.Y(frameSet);

% Assemble these components into the args variable
args = {perimeter, glintData, ellipseRMSE, gazeTargets};

% Assemble the key-values
keyVals = {...
    'eyePoseLB', p.Results.eyePoseLB,...
    'eyePoseUB', p.Results.eyePoseUB,...
    };


%% Set up the parallel pool
if p.Results.useParallel
    startParpool( p.Results.nWorkers, p.Results.verbose );
end


%% Set x0
x0 = p.Results.sceneParamsX0;
x = x0;


%% Define BADS search options
options = bads('defaults');          % Get a default OPTIONS struct
options.Display = 'off';             % Silence display output
options.UncertaintyHandling = 0;     % The objective is deterministic


%% Create a directory to save the diagnostic plots
[sceneGeomPath,sceneGeomName,] = fileparts(sceneGeometryFileName);
diagnosticDirName = fullfile(sceneGeomPath,[sceneGeomName '_diagnostics']);
if ~exist(diagnosticDirName, 'dir')
    mkdir(diagnosticDirName);
else
    rmdir(diagnosticDirName, 's');
    mkdir(diagnosticDirName);
end


%% Loop over iterations

for ii = 1:p.Results.searchIterations
    
    %% Set up the fit figure
    nStages = 4;
    figHandle = addSubPlots([],0,nStages);
    boundTol = 1e-6;
    addPlotsWrap = @(idx,x,fitAtBound) addSubPlots(figHandle,idx,nStages,x,sceneGeometry,args{:},fitAtBound,keyVals);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% STAGE 1 -- TORSION / TRANSLATION SEARCH
    % Perform an initial, iterated search, locking parameters for camera
    % distance, eye rotation, and corneal curvature.
    
    % Announce
    if p.Results.verbose
        fprintf(['Iter 0' num2str(ii) ', Stage 1...']);
    end
    % Bounds
    bound = [20, 10, 10, 0, 0, 0, 0, 0, 0];
    lb = x - bound;
    ub = x + bound;
    lbp = x - bound./2;
    ubp = x + bound./2;
    [lb,ub,lbp,ubp] = cornealCurvConstraint(sceneGeometry,lb,ub,lbp,ubp);
    % Search
    x = iterativeSearch(x,sceneGeometry,args,keyVals,lb,ub,lbp,ubp,options);
    xStages(1,:) = x;
    % Identify any params that hit a bound in the final search stage
    notLocked = lb ~= ub;
    fitAtBound = any([(abs(x(notLocked)-lb(notLocked)) < boundTol); (abs(x(notLocked)-ub(notLocked)) < boundTol)]);
    % Plot
    addPlotsWrap(1,x,fitAtBound);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% STAGE 2 -- ROTATION CENTER SEARCH
    % Search over the eye rotation center
    
    % Announce
    if p.Results.verbose
        fprintf('Stage 2...');
    end
    % Bounds
    lb = [x(1:4), 0.75, 0.75, x(7:9)];
    ub = [x(1:4), 1.25, 1.25, x(7:9)];
    lbp = [x(1:4), 0.85, 0.85, x(7:9)];
    ubp = [x(1:4), 1.15, 1.15, x(7:9)];
    [lb,ub,lbp,ubp] = cornealCurvConstraint(sceneGeometry,lb,ub,lbp,ubp);
    % Objective
    myObj = @(x) calcGlintGazeError( updateSceneGeometry( sceneGeometry, x ), args{:}, keyVals{:} );
    % Search
    x = bads(myObj,x,lb,ub,lbp,ubp,[],options);
    xStages(2,:) = x;
    % Identify any params that hit a bound in the final search stage
    notLocked = lb ~= ub;
    fitAtBound = any([(abs(x(notLocked)-lb(notLocked)) < boundTol); (abs(x(notLocked)-ub(notLocked)) < boundTol)]);
    % Plot
    addPlotsWrap(2,x,fitAtBound);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% STAGE 3 -- TRANSLATION AND CURVATURE SEARCH
    % Lock the rotation centers, search over translation and corneal curvature
    
    % Announce
    if p.Results.verbose
        fprintf('Stage 3...');
    end
    % Bounds
    bound = [abs(x(1:3).*0.25), 0, 0, 0, x(7:8).*0.25 90];
    lb = x - bound;
    ub = x + bound;
    lbp = x - bound./2;
    ubp = x + bound./2;
    [lb,ub,lbp,ubp] = cornealCurvConstraint(sceneGeometry,lb,ub,lbp,ubp);
    % Search
    x = iterativeSearch(x,sceneGeometry,args,keyVals,lb,ub,lbp,ubp,options);
    xStages(3,:) = x;
    % Identify any params that hit a bound in the final search stage
    notLocked = lb ~= ub;
    fitAtBound = any([(abs(x(notLocked)-lb(notLocked)) < boundTol); (abs(x(notLocked)-ub(notLocked)) < boundTol)]);
    % Plot
    addPlotsWrap(3,x,fitAtBound);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% STAGE 4 -- COMPLETE SEARCH
    % Search over all parameters
    
    % Announce
    if p.Results.verbose
        fprintf('Stage 4...');
    end
    % Bounds
    bb = 0.1 / ii;
    lb  = [x(1:8)./((1-bb).^-sign(x(1:8))), x(9)-10];
    lbp = [x(1:8)./((1-bb/2).^-sign(x(1:8))), x(9)-5];
    ubp = [x(1:8)./((1+bb/2).^-sign(x(1:8))), x(9)+5];
    ub  = [x(1:8)./((1+bb).^-sign(x(1:8))), x(9)+10];
    [lb,ub,lbp,ubp] = cornealCurvConstraint(sceneGeometry,lb,ub,lbp,ubp);
    % Lock the depth parameter if so instructed
    if p.Results.lockDepth
        lb(4) = x(4); lbp(4) = x(4); ubp(4) = x(4); ub(4) = x(4);
    end
    % Objective
    myObj = @(x) calcGlintGazeError( updateSceneGeometry( sceneGeometry, x ), args{:}, keyVals{:} );
    % Search
    x = bads(myObj,x,lb,ub,lbp,ubp,[],options);
    xStages(4,:) = x;
    % Identify any params that hit a bound in the final search stage
    notLocked = lb ~= ub;
    fitAtBound = any([(abs(x(notLocked)-lb(notLocked)) < boundTol); (abs(x(notLocked)-ub(notLocked)) < boundTol)]);
    % Plot
    addPlotsWrap(4,x,fitAtBound);
    
    
    %% Save the fit-by-stage plot
    if p.Results.verbose
        fprintf('Saving plots...\n');
    end
    
    % Save the staged fit results
    figureName = fullfile(diagnosticDirName,[sceneGeomName '_fitsByStage_iter0' num2str(ii) '.pdf']);
    addSupTitle(figHandle,sceneGeomName);
    saveas(figHandle,figureName)
    
    
end % searchIterations


%% Update the sceneGeometry 

% The fitted sceneGeometry
f = updateSceneGeometry( sceneGeometry, x );

% Obtain the model components at the solution
[ fVal, modelEyePose, modelPupilEllipse, modelGlintCoord, modelPoseGaze, modelVecGaze, poseRegParams, vectorRegParams, rawErrors] = ...
    calcGlintGazeError( f, args{:}, keyVals{:} );

% The varargin originally used to create the sceneGeometry
sceneGeometryVarargin = sceneGeometry.meta.createSceneGeometry.varargin;

% keys and values to update
keys = {...
    'cameraTorsion',...
    'cameraTranslation',...
    'rotationCenters',...
    'measuredCornealCurvature',...
    'fixationEyePose',...
    'screenTorsion',...
    'screenRotMat',...
    };
values = {...
    f.cameraPosition.torsion, ...
    f.cameraPosition.translation, ...
    f.eye.rotationCenters, ...
    f.eye.cornea.kvals, ...
    poseRegParams.t, ...
    poseRegParams.theta, ...    
    poseRegParams.R, ...
    };

% Loop through the keys and either update or add
for kk = 1:length(keys)
    idx = find(strcmp(sceneGeometryVarargin,keys{kk}),1);
    if isempty(idx)
        sceneGeometryVarargin = [sceneGeometryVarargin, keys{kk}, values{kk}];
    else
        sceneGeometryVarargin(idx+1)=values(kk);
    end
end

% Get the execution time
executionTime = toc(ticObject);

% Create a new sceneGeometry with the update key-values
sceneGeometry = createSceneGeometry(sceneGeometryVarargin{:});

% Update the meta data
sceneGeometry.meta.estimateSceneParams = p.Results;
sceneGeometry.meta.estimateSceneParams.x0 = x0;
for ii = 1:nStages
    sceneGeometry.meta.estimateSceneParams.(['x' num2str(ii)]) = xStages(ii,:);
end
sceneGeometry.meta.estimateSceneParams.fitAtBound = fitAtBound;
sceneGeometry.meta.estimateSceneParams.fVal = fVal;
sceneGeometry.meta.estimateSceneParams.executionTime = executionTime;
sceneGeometry.meta.estimateSceneParams.varargin = varargin;
sceneGeometry.meta.estimateSceneParams.sceneGeometryVarargin = sceneGeometryVarargin;
sceneGeometry.meta.estimateSceneParams.modelEyePose = modelEyePose;
sceneGeometry.meta.estimateSceneParams.modelPupilEllipse = modelPupilEllipse;
sceneGeometry.meta.estimateSceneParams.modelGlintCoord = modelGlintCoord;
sceneGeometry.meta.estimateSceneParams.modelPoseGaze = modelPoseGaze;
sceneGeometry.meta.estimateSceneParams.modelVecGaze = modelVecGaze;
sceneGeometry.meta.estimateSceneParams.poseRegParams = poseRegParams;
sceneGeometry.meta.estimateSceneParams.vectorRegParams = vectorRegParams;
sceneGeometry.meta.estimateSceneParams.rawErrors = rawErrors;

% Save the sceneGeometry file
if ~isempty(sceneGeometryFileName)
    save(sceneGeometryFileName,'sceneGeometry');
end


%% Save the final diagnostic model fit montage

% Find the video for this pupil file
if ~isempty(p.Results.grayVideoName)
    grayVideoName = p.Results.grayVideoName;
else
    grayVideoName = strrep(pupilFileName,p.Results.pupilFileToVideoSuffixSwitch{1},p.Results.pupilFileToVideoSuffixSwitch{2});
end

% Get the modeled eye poses
[ ~, modelEyePose] = calcGlintGazeError( sceneGeometry, args{:}, keyVals{:} );

% Create an eye model montage
figureName = fullfile(diagnosticDirName,[sceneGeomName '_sceneDiagnosticMontage_eyeModel.png']);
saveEyeModelMontage(sceneGeometry, modelEyePose, frameSet, grayVideoName, figureName);


%% alert the user that we are done with the routine
if p.Results.verbose
    fprintf([num2str(executionTime) '\n']);
end

end






%%%%%%%%%%%% LOCAL FUNCTIONS


function [lb,ub,lbp,ubp] = cornealCurvConstraint(sceneGeometry,lb,ub,lbp,ubp)
% The first of the corneal curvature values must always be smaller than the
% second. This constrains the differential scaling value that can be
% present for the last parameter.
kvals = sceneGeometry.eye.cornea.kvals;
diffScaleUpperBound = kvals(2)/kvals(1);
lb(8) = min([lb(8) diffScaleUpperBound]);
ub(8) = min([ub(8) diffScaleUpperBound]);
lbp(8) = min([lbp(8) diffScaleUpperBound]);
ubp(8) = min([ubp(8) diffScaleUpperBound]);

% Keep the corneal axis between +-90
lb(9) = max([lb(9) -90]);
ub(9) = min([ub(9) 90]);
lbp(9) = max([lbp(9) -90]);
ubp(9) = min([ubp(9) 90]);
end


function [x, fVal] = iterativeSearch(x,sceneGeometry,args,keyVals,lb,ub,lbp,ubp,options)
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
    [ ~, modelEyePose] = calcGlintGazeError( updateSceneGeometry( sceneGeometry, x ), args{:}, keyVals{:});
    % Update the objective
    myObj = @(x) calcGlintGazeError( updateSceneGeometry( sceneGeometry, x ), args{:}, keyVals{:}, 'modelEyePose', modelEyePose );
    % Perform the search
    x = bads(myObj,x,lb,ub,lbp,ubp,[],options);
    % The objective we care about is for the complete objective
    fVal = calcGlintGazeError( updateSceneGeometry( sceneGeometry, x ), args{:}, keyVals{:} );
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


function addSupTitle(figHandle,str)
    set(0, 'CurrentFigure', figHandle);
    gcf;
    axes('Position',[0 0 1 1],'visible','off','Tag','suptitle');
    ht=text(.5,0.98,str);set(ht,'horizontalalignment','center','fontsize',14,'Interpreter','none');
    drawnow
end
    
    


function figHandle = addSubPlots(figHandle,idx,nStages,x,sceneGeometry,perimeter,glintData,ellipseRMSE,gazeTargets,fitAtBound,keyVals)

% Prepare the figure
if idx == 0
    figHandle=figure('Visible','off');
    set(gcf,'PaperOrientation','landscape');
    set(figHandle, 'Units','inches')
    height = 11;
    width = 11;
    
    % The last two parameters of 'Position' define the figure size
    set(figHandle, 'Position',[25 5 width height],...
        'PaperSize',[width height],...
        'PaperPositionMode','auto',...
        'Color','w');
    return
else
    set(0, 'CurrentFigure', figHandle)
end

% Get the model output
[ ~, ~, modelPupilEllipse, modelGlintCoord, modelPoseGaze, modelVecGaze, ~, ~, rawErrors] = ...
    calcGlintGazeError( updateSceneGeometry( sceneGeometry, x ), perimeter, glintData, ellipseRMSE, gazeTargets, keyVals{:});

% We are going to have four sub-plots
nCols = 4;


% 4. Glint-pupil vec matching gaze targets
subplot(nStages,nCols,(idx-1)*nCols+4)
plot(gazeTargets(1,:),gazeTargets(2,:),'ok'); hold on;
plot(modelVecGaze(1,:),modelVecGaze(2,:),'xr'); hold on;
ylim([-10 10])
axis equal
myLabel = sprintf('Gaze vec [%2.2f]',rawErrors(4));
title(myLabel);

% 3. EyePose matching gaze targets
subplot(nStages,nCols,(idx-1)*nCols+3)
plot(gazeTargets(1,:),gazeTargets(2,:),'ok'); hold on;
plot(modelPoseGaze(1,:),modelPoseGaze(2,:),'xr'); hold on;
ylim([-10 10])
axis equal
myLabel = sprintf('Gaze pose [%2.2f]',rawErrors(3));
title(myLabel);

% 2. Glint fits
subplot(nStages,nCols,(idx-1)*nCols+2)
plot(glintData.X,glintData.Y,'ok'); hold on;
plot(modelGlintCoord.X,modelGlintCoord.Y,'xr'); hold on;
axis equal
myLabel = sprintf('Glint [%2.2f]',rawErrors(2));
title(myLabel);

% 1. Perimeter fits
% Define a figure
hFig = figure( 'Visible', 'off');
dim = 150;
imshow(ones(dim,dim),'Border', 'tight');
drawnow;
hAxes = get(hFig,'CurrentAxes');
hold on;
for ii = 1:length(ellipseRMSE)
    Xp = perimeter.data{ii}.Xp;
    meanXp = mean(Xp);
    Xp = Xp - meanXp + dim/2;
    Yp = perimeter.data{ii}.Yp;
    meanYp = mean(Yp);
    Yp = Yp - meanYp + dim/2;
    p1 = plot(hAxes,Xp,Yp,'.k');
    xlim([1 dim]);
    ylim([1 dim]);
    drawnow;
    hold on;
    pupilEllipseParams = modelPupilEllipse(ii,:);
    pupilEllipseParams(1) = pupilEllipseParams(1) - meanXp + dim/2;
    pupilEllipseParams(2) = pupilEllipseParams(2) - meanYp + dim/2;
    p2 = addTransparentEllipseToFigure(pupilEllipseParams,dim,dim,'red',1,hAxes);
    axis off;
    drawnow;
    thisFrame=getframe(hFig);
    framesToMontage(:,:,:,ii) = thisFrame.cdata;
    delete(p1); delete(p2);
    drawnow;
end
close(hFig)
set(0, 'CurrentFigure', figHandle)
subplot(nStages,nCols,(idx-1)*nCols+1)
montage(framesToMontage)
myLabel = sprintf('Perimeter [%2.2f]',rawErrors(1));
title(myLabel);

% Text label that indicates stage
xRange=get(gca,'XLim');
yRange=get(gca,'YLim');
ht = text(0*xRange(1)-0.2*xRange(2),0.5*yRange(2),['Stage ' num2str(idx)]);
set(ht,'Rotation',90)
set(ht,'FontSize',18)
drawnow

% If this is the last panel, put an annotation for the x parameters at the
% bottom. Report params that hit a bound in red.
if idx == nStages
    gcf;
    axes('Position',[0 0 1 1],'Visible','off','Tag','subtitle');
    str = sprintf('Camera torsion: $color-start$%2.1f$$color-end$$, position: [$color-start$%2.1f$$color-end$$, $color-start$%2.1f$$color-end$$, $color-start$%2.1f$$color-end$$]; Rotation center joint, diff [$color-start$%2.2f$$color-end$$, $color-start$%2.2f$$color-end$$]; Corneal curvature joint, diff, angle [$color-start$%2.2f$$color-end$$, $color-start$%2.2f$$color-end$$, $color-start$%2.2f$$color-end$$]',x);
    tagIdx = strfind(str,'$color-start$');
    for ii=1:length(fitAtBound)
        if fitAtBound(ii)
            str(tagIdx(ii):tagIdx(ii)+12) = '\color{red$$}';
        else
            str(tagIdx(ii):tagIdx(ii)+12) = '\color{black}';
        end
    end
    str = strrep(str,'$$color-end$$','\color{black}');
    str = strrep(str,'$$','');
    ht=text(.5,0.055,str);
    set(ht,'horizontalalignment','center','fontsize',12);
    str = sprintf('x = [ %2.3f, %2.3f, %2.3f, %2.3f, %2.3f, %2.3f, %2.3f, %2.3f, %2.3f ]',x);
    ht=text(.5,0.025,str);set(ht,'horizontalalignment','center','fontsize',12);
    drawnow
end

end




function saveEyeModelMontage(sceneGeometry, modelEyePose, frameSet, grayVideoName, montageFileName)
% Saves a montage with the model eye superimposed.

% Silence some errors that can arise during the forward projection
warningState = warning;
warning('off','projectModelEye:ellipseFitFailed');

% Sort the ellipse array list so that the frames appear in temporal order
[frameSet, sortOrder] = sort(frameSet);
modelEyePose = modelEyePose(sortOrder,:);

% Check that the file exists
if exist(grayVideoName,'file') && ~isempty(frameSet)
    
    % Open the video object
    videoInObj = videoIOWrapper(grayVideoName,'ioAction','read');
    
    % Get the video properties
    videoSizeX = videoInObj.Width;
    videoSizeY = videoInObj.Height;
    
    % Define a variable to hold the selected frames
    framesToMontage = zeros(videoSizeY,videoSizeX,3,length(frameSet),'uint8');
    
    % Define a figure
    hFig = figure( 'Visible', 'off');
    hAxes = gca();
    
    % Loop through the frames and keep the matching ones
    for ii = 1:length(frameSet)
        idx = frameSet(ii);
        videoInObj.CurrentTime = (idx - 1)/(videoInObj.FrameRate);
        sourceFrame = readFrame(videoInObj);
        imshow(sourceFrame,'Border', 'tight','Parent',hAxes);
        hold on
        axis off;
        % Add the rendered eye model
        eyePose = modelEyePose(ii,:);
        if ~any(isnan(eyePose))
            renderEyePose(eyePose, sceneGeometry, 'newFigure', false, ...
                'modelEyeLabelNames', {'retina' 'irisPerimeter' 'pupilEllipse' 'cornea' 'glint_01' 'glint_02'}, ...
                'modelEyePlotColors', {'.w' '.b' '-g' '.y' 'xr' 'xr'}, ...
                'modelEyeAlpha', [0.25 0.25 0.25 0.25 1 1],...
                'modelEyeSymbolSizeScaler',1.5,...
                'showAzimuthPlane',true);
        end
        % Get the frame
        drawnow;
        thisFrame=getframe(hFig);
        % Add a text label for the frame number
        frameLabel = sprintf('frame: %d',idx);
        thisFrame.cdata = insertText(thisFrame.cdata,[20 20],frameLabel,'FontSize',30);
        % Store the frame. Detect if we have a bad or empty frame and then
        % skip if that is the case
        if all(size(squeeze(framesToMontage(:,:,:,ii)))==size(thisFrame.cdata))
            framesToMontage(:,:,:,ii) = thisFrame.cdata;
        end
        % hold off
        hold off
    end
    
    % Close the temporary figure
    close(hFig);
    
    % Prepare the figure
    figHandle=figure('Visible','off');
    set(gcf,'PaperOrientation','landscape');
    set(figHandle, 'Units','inches')
    height = 6;
    width = 11;
    
    % The last two parameters of 'Position' define the figure size
    set(figHandle, 'Position',[25 5 width height],...
        'PaperSize',[width height],...
        'PaperPositionMode','auto',...
        'Color','w');
    
    % Turn off a warning that can occur during the montage step
    warningState = warning;
    warning('off','images:initSize:adjustingMag');
    
    % Create the montage
    montage(framesToMontage);
    
    % Restore the warning state
    warning(warningState);
    
    % Save the montage
    saveas(figHandle,montageFileName)
    
    % Close the figure
    close(figHandle)
    
    % Rotate the figure by 90 degrees clockwise, because I can't get the
    % MATLAB plotting routines to output the image how I want it.
    A = imread(montageFileName);
    A = rot90(A,3);
    imwrite(A,montageFileName);
    
    % close the video object
    clear videoInObj
    
end % There is a file to plot

% Restore the warning state
warning(warningState);

end % saveEyeModelMontage




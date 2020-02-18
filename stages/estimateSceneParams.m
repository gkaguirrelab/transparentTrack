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
%   also necessary to provide a reasonably accurate value for the distance
%   of the camera from the eye.
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
    estimateSceneParams(pupilFileName, perimeterFileName, glintFileName, sceneGeometryFileName, 'frameSet', frameSet, 'gazeTargets', gazeTargets);
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
p.addParameter('sceneParamsX0',[0 0 0 130 1 1 1 1],@isnumeric);
p.addParameter('lockDepth',true,@islogical);
p.addParameter('eyePoseLB',[-89,-89,0,0.1],@isnumeric);
p.addParameter('eyePoseUB',[89,89,0,4],@isnumeric);
p.addParameter('frameSet',[],@(x)(isempty(x) | isnumeric(x)));
p.addParameter('gazeTargets',[],@(x)(isempty(x) | isnumeric(x)));
p.addParameter('fixSpectacleLens',[],@(x)(isempty(x) | isnumeric(x)));

% parse
p.parse(pupilFileName, perimeterFileName, glintFileName, sceneGeometryFileName, varargin{:})




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


%% Set up the fit figure
nStages = 4;
figHandle = addSubPlots([],0,nStages);
addPlotsWrap = @(idx,x) addSubPlots(figHandle,idx,nStages,x,sceneGeometry,args{:},keyVals);


%% Set x0
x0 = p.Results.sceneParamsX0;
x = x0;

%% Define BADS search options
options = bads('defaults');          % Get a default OPTIONS struct
options.Display = 'off';             % Silence display output
options.UncertaintyHandling = 0;     % The objective is deterministic


%% STAGE 1 -- TORSION / TRANSLATION SEARCH
% Perform an initial, iterated search, locking parameters for camera
% distance, eye rotation, and corneal curvature.

% Announce
if p.Results.verbose
    fprintf('Stage 1...');
end
% Bounds
bound = [20, 10, 10, 0, 0, 0, 0, 0];
lb = x - bound;
ub = x + bound;
lbp = x - bound./2;
ubp = x + bound./2;
% Search
x = iterativeSearch(x,sceneGeometry,args,keyVals,lb,ub,lbp,ubp,options);
xStages(1,:) = x;
% Plot
addPlotsWrap(1,x);


%% STAGE 2 -- ROTATION CENTER SEARCH
% Search over the eye rotation center

% Announce
if p.Results.verbose
    fprintf('Stage 2...');
end
% Bounds
lb = [x(1:4), 0.75, 0.75, x(7:8)];
ub = [x(1:4), 1.25, 1.25, x(7:8)];
lbp = [x(1:4), 0.85, 0.85, x(7:8)];
ubp = [x(1:4), 1.15, 1.15, x(7:8)];
% Objective
myObj = @(x) calcGlintGazeError( updateSceneGeometry( sceneGeometry, x ), args{:}, keyVals{:} );
% Search
x = bads(myObj,x,lb,ub,lbp,ubp,[],options);
xStages(2,:) = x;
% Plot
addPlotsWrap(2,x);


%% STAGE 3 -- TRANSLATION AND CURVATURE SEARCH
% Lock the rotation centers, search over translation and corneal curvature

% Announce
if p.Results.verbose
    fprintf('Stage 3...');
end
% Bounds
bound = [abs(x(1:3).*0.25), 0, 0, 0, x(7:8).*0.25];
lb = x - bound;
ub = x + bound;
lbp = x - bound./2;
ubp = x + bound./2;
% Search
x = iterativeSearch(x,sceneGeometry,args,keyVals,lb,ub,lbp,ubp,options);
xStages(3,:) = x;
% Plot
addPlotsWrap(3,x);


%% STAGE 4 -- COMPLETE SEARCH
% Search over all parameters

% Announce
if p.Results.verbose
    fprintf('Stage 4...');
end
% Bounds
lb  = x./(0.90.^-sign(x));
lbp = x./(0.95.^-sign(x));
ubp = x./(1.05.^-sign(x));
ub  = x./(1.10.^-sign(x));
% if we have been told to lock the depth parameter, do so
if p.Results.lockDepth
    lb(4) = x(4); lbp(4) = x(4); ubp(4) = x(4); ub(4) = x(4);
end
% Objective
myObj = @(x) calcGlintGazeError( updateSceneGeometry( sceneGeometry, x ), args{:}, keyVals{:} );
% Search
x = bads(myObj,x,lb,ub,lbp,ubp,[],options);
xStages(4,:) = x;
% Plot
addPlotsWrap(4,x);


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

% Update and move the meta data around
sceneGeometry.meta.estimateSceneParams = p.Results;
sceneGeometry.meta.estimateSceneParams.x0 = x0;
for ii = 1:nStages
    sceneGeometry.meta.estimateSceneParams.(['x' num2str(ii)]) = xStages(ii,:);
end
sceneGeometry.meta.estimateSceneParams.fVal = fVal;
sceneGeometry.meta.estimateSceneParams.executionTime = executionTime;
sceneGeometry.meta.estimateSceneParams.varargin = varargin;
sceneGeometry.meta.estimateSceneParams.sceneGeometryVarargin = sceneGeometryVarargin;

% Add the model components at the solution to the meta data
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


%% Save diagnostic plots

if p.Results.verbose
    fprintf('Saving plots.\n');
end
[sceneGeomPath,sceneGeomName,] = fileparts(sceneGeometryFileName);
diagnosticDirName = fullfile(sceneGeomPath,[sceneGeomName '_diagnostics']);
if ~exist(diagnosticDirName, 'dir')
    mkdir(diagnosticDirName);
else
    rmdir(diagnosticDirName, 's');
    mkdir(diagnosticDirName);
end

% Find the video for this pupil file
if ~isempty(p.Results.grayVideoName)
    grayVideoName = p.Results.grayVideoName;
else
    grayVideoName = strrep(pupilFileName,p.Results.pupilFileToVideoSuffixSwitch{1},p.Results.pupilFileToVideoSuffixSwitch{2});
end

% Get the modeled eye poses
[ ~, modelEyePose] = calcGlintGazeError( sceneGeometry, args{:}, keyVals{:} );

% Save the staged fit results
figureName = fullfile(diagnosticDirName,[sceneGeomName '_fitsByStage.pdf']);
addSupTitle(figHandle,sceneGeomName);
saveas(figHandle,figureName)

% Create an eye model montage
figureName = fullfile(diagnosticDirName,[sceneGeomName '_sceneDiagnosticMontage_eyeModel.png']);
saveEyeModelMontage(sceneGeometry, modelEyePose, frameSet, grayVideoName, figureName)


%% alert the user that we are done with the routine
if p.Results.verbose
    executionTime
    fprintf('\n');
end

end






%%%%%%%%%%%% LOCAL FUNCTIONS




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
    gcf
    axes('Position',[0 0 1 1],'visible','off','Tag','suptitle');
    ht=text(.5,0.98,str);set(ht,'horizontalalignment','center','fontsize',14,'Interpreter','none');
    drawnow
end
    
    


function figHandle = addSubPlots(figHandle,idx,nStages,x,sceneGeometry,perimeter,glintData,ellipseRMSE,gazeTargets,keyVals)

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
% bottom.
if idx == nStages
    gcf
    axes('Position',[0 0 1 1],'visible','off','Tag','subtitle');
    str = sprintf('Camera torsion: %2.1f, position: [%2.1f, %2.1f, %2.1f]; Rotation center joint, diff [%2.2f, %2.2f]; Corneal curvature joint diff [%2.2f, %2.2f]',x);
    ht=text(.5,0.05,str);set(ht,'horizontalalignment','center','fontsize',12);
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
    figHandle=figure('visible','off');
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




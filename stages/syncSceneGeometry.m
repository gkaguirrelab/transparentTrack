function syncSceneGeometry(pupilFileName, varargin)
% Update sceneGeometry camera position to match an acquisition
%
% Syntax:
%  syncSceneGeometry(pupilFileName, varargin)
%
% Description:
%   A sceneGeometry file is created for a given acquisition. Included in
%   the sceneGeometry is a specification of properties of the extrinsic
%   camera matrix, including the position of the camera in space relative
%   to the coordinates, which have as their origin the anterior surface of
%   the cornea along the optical axis of the eye. If we wish to use this
%   sceneGeometry file for the analysis of data from other acqusitions for
%   a given subject, we need to deal with the possibility that the subject
%   has moved their head between acquisitions. As the coordinate system is
%   based upon a fixed anatomical landmark of the eye, the effect of head
%   translation in this system is to change the camera position. This
%   routine assists in calculating an updated camera position for a given
%   acquisition.
%
%   Select a sceneGeometry with the UI file picker, then select the
%   acquisitions to which you would like to align the sceneGeometry.
%   An image derived from the sceneGeometry and the acquisition is shown,
%   and controls are used to adjust the acquisition to match the
%   sceneGeometry.
%
% Inputs:
%   pupilFileName         - Full path to the .mat file that contains the
%                           pupil data to which the sceneGeometry should be
%                           synced.
%   
% Optional key/value pairs (display and I/O):
%  'sceneGeometryFileNameToSync' - Full path to the .mat file that contains
%                           the sceneGeometry to be used.
%  'verbose'              - Logical. Default false.
%  'displayMode'          - Logical. Default false.
%  'saveAdjustedSceneGeometry' - Logical.
%  'saveDiagnosticPlot'   - Logical.
%  'doNotSyncSceneToItself' - Logical.
%   
% Optional key/value pairs (fitting params):
%  'alignMethod'          - Char vec, or cell array of char vecs.
%  'deltaPix'             - Numeric.
%  'deltaDeg'             - Numeric.
%  'deltaScale'           - Numeric.
%  'deltaPose'            - Numeric.
%  'eyePositionTargetLengthFrames' - Scalar.
%  'gazeErrorThreshTol'   - Scalar.
%
% Examples:
%{
    % Invoke the file picker GUI
    syncSceneGeometry('','displayMode',true,'verbose',true,'alignMethod','shape');
%}
%{
    sceneGeometryFileNameToSync = '/Users/aguirre/Dropbox (Aguirre-Brainard Lab)/TOME_processing/session2_spatialStimuli/TOME_3015/032417/EyeTracking/GazeCal03_sceneGeometry.mat';
    pupilFileName = '/Users/aguirre/Dropbox (Aguirre-Brainard Lab)/TOME_processing/session2_spatialStimuli/TOME_3015/032417/EyeTracking/tfMRI_MOVIE_AP_run01_pupil.mat';
    syncSceneGeometry(pupilFileName, 'sceneGeometryFileNameToSync', sceneGeometryFileNameToSync) 
%}

%% Parse vargin for options passed here
p = inputParser; p.KeepUnmatched = true;

% Required
p.addRequired('pupilFileName',@ischar);

% Optional display and I/O params
p.addParameter('sceneGeometryFileNameToSync','',@(x)(ischar(x) | iscell(x)));
p.addParameter('verbose',false,@islogical);
p.addParameter('displayMode',false,@islogical);
p.addParameter('saveAdjustedSceneGeometry',true,@islogical);
p.addParameter('saveDiagnosticPlot',true,@islogical);
p.addParameter('doNotSyncSceneToItself',true,@islogical);

% Optional flow control params
p.addParameter('useParallel',true,@islogical);
p.addParameter('nWorkers',[],@(x)(isempty(x) || isnumeric(x)));

% Optional fitting params
p.addParameter('alignMethod','gazePre',@(x)(ischar(x) | iscell(x)));
p.addParameter('sceneSyncBound',[0, 10, 10, 0, 0, 0, 0, 0],@isnumeric);
p.addParameter('sceneSyncX',[],@isnumeric);
p.addParameter('eyePositionTargetLengthFrames',10,@isscalar);
p.addParameter('gazeErrorThreshTol',0.25,@isscalar);


%% Parse and check the parameters
p.parse(pupilFileName, varargin{:});


%% Load sceneGeometryIn files
% This sceneGeometry--and its associated video, perimeter, and pupil
% data--constitute the "moving" measurements.

% Get the name and path
if iscell(p.Results.sceneGeometryFileNameToSync)
    sceneGeometryFileNameToSync = p.Results.sceneGeometryFileNameToSync{1};
else
    sceneGeometryFileNameToSync = p.Results.sceneGeometryFileNameToSync;
end
if isempty(sceneGeometryFileNameToSync)
    % Open a file picker UI to select a sceneGeometry
    [sceneGeometryInName,sceneGeometryInPath] = uigetfile(fullfile('.','*_sceneGeometry.mat'),'Choose a sceneGeometry file');
    sceneGeometryFileNameToSync = fullfile(sceneGeometryInPath,sceneGeometryInName);
else
    tmp = strsplit(sceneGeometryFileNameToSync,filesep);
    sceneGeometryInName = tmp{end};
    % Handle the case in which only the name, and not the path, to the
    % sceneGeometryIn file is specified. In this case, assume that the path
    % is the same as that given for the pupilFileName
    if length(tmp)==1
        tmp = strsplit(pupilFileName,filesep);
        sceneGeometryInPath = fullfile(filesep,tmp{1:end-1});
        sceneGeometryFileNameToSync = fullfile(sceneGeometryInPath,sceneGeometryInName);
    else
        sceneGeometryInPath = fullfile(filesep,tmp{1:end-1});
    end
end

% Load the sceneGeometry variable into memory
dataLoad=load(sceneGeometryFileNameToSync);
sceneGeometryIn=dataLoad.sceneGeometry;
clear dataLoad

% Get the filename stem for this sceneGeometry file, so that we can use it
% to load other, associated files
tmp = strsplit(sceneGeometryInName,'_sceneGeometry.mat');
sceneGeometryInStem = tmp{1};

% Load the pupil data file associated with the sceneGeometry
tmp = fullfile(sceneGeometryInPath,[sceneGeometryInStem,'_pupil.mat']);
load(tmp,'pupilData');

% Check that we have a radiusSmoothed field
if ~isfield(pupilData,'radiusSmoothed')
    warning('This sceneGeometry pupilData file does not contain a radiusSmoothed field; returning')
    return
end

% Load the perimeter file associated with the sceneGeometry
tmp = fullfile(sceneGeometryInPath,[sceneGeometryInStem '_correctedPerimeter.mat']);
load(tmp,'perimeter');


%% Derive properties from sceneGeometryIn
% The sceneGeometryIn specifies a particular eyePose as corresponding to a
% screen fixation position of [0 0]. We identify pupil, perimeter, and
% video frame measurements at this position.

% Get the camera offset point. This is where the center of the video image
% lands on the camera sensor array.
cameraOffsetPoint = [sceneGeometryIn.cameraIntrinsic.matrix(1,3), ...
    sceneGeometryIn.cameraIntrinsic.matrix(2,3)];


%% Find the fixation frame for the MOVING image
% This is the sceneGeometry file that (typically) was derived during a
% gazeCalibration procedure.

% Which of the list of frames is the [0,0] fixation frame
idx = find((sceneGeometryIn.meta.estimateSceneParams.gazeTargets(1,:)==0).*(sceneGeometryIn.meta.estimateSceneParams.gazeTargets(2,:)==0));
referenceFrameMoving = sceneGeometryIn.meta.estimateSceneParams.frameSet(idx);

% Obtain and store the pupil perimeter points for this frame
XpMoving = perimeter.data{referenceFrameMoving}.Xp;
YpMoving = perimeter.data{referenceFrameMoving}.Yp;

% Store the pupilEllipse for this frame
pupilEllipseMoving = sceneGeometryIn.meta.estimateSceneParams.modelPupilEllipse(idx,:);

% The eyePose for this reference frame is the modeled value
eyePoseMoving = sceneGeometryIn.meta.estimateSceneParams.modelEyePose(idx,:);

% Load in the image for this frame. This is the "moving" frame.
tmp = fullfile(sceneGeometryInPath,[sceneGeometryInStem '_gray.avi']);
videoFrameMoving = makeMedianVideoImage(tmp,'startFrame',referenceFrameMoving,'nFrames',1);

% Find the shape of the pupil for this frame, expressed as theta and rho
% values (SEE: csaEllipseError)
pupilRhoShapeFixed = 1-sqrt(1-pupilEllipseMoving(4)^2);
pupilThetaShapeFixed = pupilEllipseMoving(5)*2;


%% Load sceneGeometryOut files
% This acquisition--and its associated video, perimeter, and pupil
% data--constitute the "fixed" measurements. Our goal is to adjust the
% sceneGeometryIn to best match this acquisition.

% If the pupilFileName is not defined, offer some choices
if isempty(pupilFileName)
    % Get a list of all gray.avi videos in this directory
    fileList = dir(fullfile(sceneGeometryInPath,'*_pupil.mat'));
    
    % Exclude the video that is the source of the fixed image
    keep=cellfun(@(x) ~strcmp(x,[sceneGeometryInStem '_pupil.mat']),extractfield(fileList,'name'));
    fileList = fileList(keep);
    
    % Ask the operator which of the videos we wish to adjust
    fprintf('\n\nSelect the pupil data to target:\n')
    for pp=1:length(fileList)
        optionName=['\t' num2str(pp) '. ' fileList(pp).name '\n'];
        fprintf(optionName);
    end
    fprintf('\nEnter a single acquisition number:\n')
    choice = input('\nYour choice: ','s');
    fileList = fileList(eval(choice));
    
    tmp = strsplit(fileList(1).name,'_pupil.mat');
    sceneGeometryOutStem = tmp{1};
    sceneGeometryOutPath = fileList(1).folder;
else
    % Derive the path and file name stem for this acquisition
    tmp = strsplit(pupilFileName,filesep);
    sceneGeometryOutPath = fullfile(filesep,tmp{1:end-1});
    tmp = strsplit(tmp{end},'_pupil.mat');
    sceneGeometryOutStem = tmp{1};
end

% Load the timebase, pupilData, perimeter, and relative camera position for
% this acquisition
tmp = fullfile(sceneGeometryOutPath,[sceneGeometryOutStem '_timebase.mat']);
load(tmp,'timebase');
tmp =  fullfile(sceneGeometryOutPath,[sceneGeometryOutStem '_pupil.mat']);
load(tmp,'pupilData');
tmp =  fullfile(sceneGeometryOutPath,[sceneGeometryOutStem '_glint.mat']);
load(tmp,'glintData');
tmp =  fullfile(sceneGeometryOutPath,[sceneGeometryOutStem '_correctedPerimeter.mat']);
load(tmp,'perimeter');
tmp =  fullfile(sceneGeometryOutPath,[sceneGeometryOutStem '_relativeCameraPosition.mat']);
load(tmp,'relativeCameraPosition');

% Identify the acqStartTimeFixed, which is the time point at which the
% fMRI acquisition began
[~, acqStartFrameFixed] = min(abs(timebase.values));


%% Check if we are syncing a sceneGeometry to itself
if p.Results.doNotSyncSceneToItself && strcmp(sceneGeometryInStem,sceneGeometryOutStem)
    if p.Results.verbose
        fprintf('Detected that the sceneGeometry source is the same as the target acquisition; returning.\n');
    end
    return
end


%% Find the target window for the acqusition
% Identify a window of frames from the acquisition during which the eye is
% fixating "screen position" [0, 0], and an error metric for the
% difference of fixation position between the sceneGeometryIn and each
% frame of the acquisition. The approach varies based upon the alignMethod
% flag.
alignMethod = p.Results.alignMethod;
if iscell(alignMethod)
    alignMethod = alignMethod{1};
end
switch alignMethod
    case 'gazePre'
        % For some acquisitions, the subject was asked to stare at a
        % fixation point in the center of the screen prior to the start of
        % the acquisition. Find the period prior to the start of the scan
        % when the eye was in the most consistent position, and closest to
        % the median position.
        windowStart = 1;
        windowEnd = acqStartFrameFixed;
        gazeX = pupilData.initial.ellipses.values(windowStart:windowEnd,1);
        gazeY = pupilData.initial.ellipses.values(windowStart:windowEnd,2);
        medianX = nanmedian(gazeX);
        medianY = nanmedian(gazeY);
        
        % Assuming that the eye had a central tendency of fixation upon the
        % center of the screen, this vector expresses the deviation of
        % fixation on any given frame from the screen center.
        eyeMatchError = sqrt(sum([gazeX-medianX; gazeY-medianY].^2,2));
        x0 = 0.5;
        
    case 'gazePost'
        % For some acqusitions (i.e., retinotopy), the subject was asked to
        % stare at a fixation point in the center of the screen after the
        % start of the acquisition. This could also work for movie viewing,
        % if we are willing to assume that the median gaze position during
        % the first 10 seconds of a movie is the center of the screen.
        % Find the period after to the start of the scan when the eye was
        % in the most consistent position, and closest to the median
        % position.
        windowStart = acqStartFrameFixed;
        windowEnd = acqStartFrameFixed+600;
        gazeX = pupilData.initial.ellipses.values(windowStart:windowEnd,1);
        gazeY = pupilData.initial.ellipses.values(windowStart:windowEnd,2);
        medianX = nanmedian(gazeX);
        medianY = nanmedian(gazeY);
        
        % Assuming that the eye had a central tendency of fixation upon the
        % center of the screen, this vector expresses the deviation of
        % fixation on any given frame from the screen center.
        eyeMatchError = sqrt(sum([gazeX-medianX; gazeY-medianY].^2,2));
        x0 = 0.5;
        
    case 'shape'
        % Find the period after the start of the scan when the pupil has a
        % shape most similar to the shape from the sceneGeometry file for
        % gaze [0 0].
        windowStart = acqStartFrameFixed;
        windowEnd = size(pupilData.initial.ellipses.values,1);
        pupilRhoShapeMoving = pupilData.initial.ellipses.values(windowStart:windowEnd,4);
        pupilRhoShapeMoving = 1-sqrt(1-pupilRhoShapeMoving.^2);
        pupilThetaShapeMoving = pupilData.initial.ellipses.values(windowStart:windowEnd,5);
        pupilThetaShapeMoving = pupilThetaShapeMoving.*2;
        
        % This vector expresses the difference in pupil shape between the
        % reference period of the sceneGeometryIn and the frames of the
        % acquisition.
        shapeDifference = ...
            sqrt(pupilRhoShapeFixed^2 + pupilRhoShapeMoving.^2 - 2*pupilRhoShapeFixed.*pupilRhoShapeMoving.*cos(pupilThetaShapeFixed-pupilThetaShapeMoving))./2;
        
        % This error vector will be minimized for frames on which the shape
        % of the pupil is most similar to the reference period of the
        % sceneGeometryIn.
        eyeMatchError = sqrt(sum(shapeDifference.^2,2));
        x0 = 0.01;
        
end

% Find the minimum fixation error threshold that results in a run of
% consecutive frames at fixation of the target length.
targetLength = p.Results.eyePositionTargetLengthFrames;

% Anonymous function to grab the indicies of when runs of frames begin
pullStartIndices = @(vec) vec(1:2:end-1);

% Anonymous function to grab the length of each run
pullRunLengths = @(vec) vec(2:2:end)-pullStartIndices(vec);

% Anonynous function that provides the lengths of runs of frames for which
% the eyeMatch error is below a threshold.
runStarts = @(thresh) find(diff([0,(eyeMatchError < thresh)',0]==1));

% Set the fzero search options
options = optimset('fzero');
options.Display = 'off';

% Adjust the targetLength as needed to achieve a threshVal below threshTol.
stillWorking = true;
while stillWorking
    % An objective function that expresses the difference of the longest
    % run length from the target run length (e.g., 30 frames long). The
    % business with the min([1e6 ...]) is to handle the case when the run
    % set is empty, and thus would otherwise return an empty variable for
    % the objective.
    myObj = @(thresh) min([1e6, (targetLength - max(pullRunLengths(runStarts(thresh))))]);
    
    % Perform the search
    threshValFixed = fzero(myObj,x0,options);
    
    % Check to see if the search has met our threshold criteria, or if we
    % have run the targetLength down as short as it can go.
    if threshValFixed < p.Results.gazeErrorThreshTol || targetLength == 1
        stillWorking = false;
    else
        targetLength = targetLength-1;
    end
end

% Check if we found a solution
if ~isfinite(threshValFixed)
    warning('Unable to find a suitable set of frames from the acquisition')
    return
end

% Find the start point of this run of frames
runLengths = pullRunLengths(runStarts(threshValFixed));
runIndices = pullStartIndices(runStarts(threshValFixed))+windowStart-1;
runLengthFixed = targetLength-myObj(threshValFixed);
startIndexFixed = runIndices(runLengths == runLengthFixed);
startIndexFixed = startIndexFixed(1);
frameSetFixed = startIndexFixed:startIndexFixed+runLengthFixed-1;

% Find the frame with the lowest ellipse RMSE during the target window
ellipseRMSE = pupilData.initial.ellipses.RMSE(frameSetFixed);
referenceFrameFixed = startIndexFixed + find(ellipseRMSE == min(ellipseRMSE)) - 1;

% Extract the frames we want
XpFixed = perimeter.data{referenceFrameFixed}.Xp;
YpFixed = perimeter.data{referenceFrameFixed}.Yp;
perimeter.data = perimeter.data(frameSetFixed);
ellipseRMSE = pupilData.initial.ellipses.RMSE(frameSetFixed);
glintCoordFixed = [glintData.X(referenceFrameFixed), glintData.Y(referenceFrameFixed)];
glintData.X = glintData.X(frameSetFixed); glintData.Y = glintData.Y(frameSetFixed);

% Load in the best image from the target period for acquisition. This is
% the "fixed" frame.
tmp = fullfile(sceneGeometryOutPath,[sceneGeometryOutStem '_gray.avi']);
videoFrameFixed = makeMedianVideoImage(tmp,'startFrame',referenceFrameFixed,'nFrames',1);


%% Create the moving sceneGeometry
% We will search across camera torsion and translation parameters to fit
% the fixed frame

% For the eyePose to be equal to the fixationPose
poseBound = [2 2 0];
eyePoseLB = [eyePoseMoving(1:3)-poseBound 0.1];
eyePoseUB = [eyePoseMoving(1:3)+poseBound 4];

% Assemble these components into the args variable
args = {perimeter, glintData, ellipseRMSE, []};

% Assemble the key-values
keyVals = {...
    'eyePoseLB', eyePoseLB,...
    'eyePoseUB', eyePoseUB,...
    };

%% Set up the parallel pool
if p.Results.useParallel
    startParpool( p.Results.nWorkers, p.Results.verbose );
end

%% Define BADS search options
options = bads('defaults');          % Get a default OPTIONS struct
options.Display = 'off';             % Silence display output
options.UncertaintyHandling = 0;     % The objective is deterministic

% Bounds
x0 = [sceneGeometryIn.meta.estimateSceneParams.x4(1:4) 1 1 1 1];
bound = p.Results.sceneSyncBound;
lb = x0 - bound;
ub = x0 + bound;
lbp = x0 - bound./2;
ubp = x0 + bound./2;
% Objective
myObj = @(x) calcGlintGazeError( updateSceneGeometry( sceneGeometryIn, x ), args{:}, keyVals{:} );
% Search
if isempty(p.Results.sceneSyncX)
    x = bads(myObj,x0,lb,ub,lbp,ubp,[],options);
else
    x = p.Results.sceneSyncX;
end


%% Optional display mode
% If we are in display mode, offer an interface for the user to manually
% adjust the delta variables

if p.Results.displayMode
        
    % Provide some instructions for the operator
    fprintf([sceneGeometryInPath '\n']);
    fprintf('Adjust horizontal and vertical camera translation with the arrow keys.\n');
    fprintf('Adjust depth camera translation with + and -.\n');
    fprintf('Adjust camera torsion with j and k.\n');
    fprintf('Switch between moving and fixed image by pressing f.\n');
    fprintf('Turn on and off perimeter display with p.\n');
    fprintf('Turn on and off model display with m.\n');
    fprintf('Press esc to exit.\n\n');
    
    % Create a figure
    figHandle = figure();
    imshow(videoFrameFixed,[],'Border','tight');
    ax = gca;
    ax.Toolbar = [];
    hold on
    text(20,30,'FIXED', 'Color', 'g','Fontsize',16);
    
    % Prepare for the loop
    showMoving = false;
    showPerimeter=false;
    showModel=false;
    stillWorking = true;
    
    % Enter the while stillWorking loop
    while stillWorking
        
        % Prepare to update the image
        hold off
        
        if showMoving
            % Work with the moving frame
            displayImage = videoFrameMoving;
            
            % Update the perimeter points
            XpDisplay = XpMoving; YpDisplay = YpMoving;
                        
        else
            % Work with the fixed frame
            displayImage = videoFrameFixed;
            
            % Update the perimeter
            XpDisplay = XpFixed; YpDisplay = YpFixed;
            
        end
        
        % Display the perimeter points
        if showPerimeter
            idx = sub2ind(size(displayImage),round(YpDisplay),round(XpDisplay));
            displayImage(idx)=255;
        end
        
        % Add the eye model
        if showModel
            % Let the user know this will take a few seconds
            text_str = 'Updating model...';
            annotHandle = addAnnotation(text_str);
            % Set the eye pose, depending upon if we are looking at the
            % fixed or moving image
            if showMoving
                sceneGeometryDisplay = sceneGeometryIn;
                eyePoseDisplay = eyePoseMoving;
            else
                sceneGeometryDisplay = updateSceneGeometry( sceneGeometryIn, x );
                eyePoseDisplay = eyePoseEllipseFit(XpFixed, YpFixed, sceneGeometryDisplay,'glintCoord',glintCoordFixed,keyVals{:});
            end
            % Render the eye model
            [tmpHandle,~,displayImage]=renderEyePose(eyePoseDisplay, sceneGeometryDisplay, ...
                'newFigure', true, 'visible', false, ...
                'backgroundImage', displayImage, ...
                'showAzimuthPlane', true, ...
                'modelEyeLabelNames', {'retina' 'pupilEllipse' 'cornea' 'glint_01'}, ...
                'modelEyePlotColors', {'.w' '-g' '.y' 'xr'}, ...
                'modelEyeSymbolSizeScaler',1.5,...
                'modelEyeAlpha', [0.25 0.25 0.25 1]);
            close(tmpHandle);
            displayImage = displayImage.cdata;
            % Remove the updating annotation
            delete(annotHandle);
        end

        % Move the moving image
        if showMoving
            % Get the 2D registration params that brings the movingImage
            % into the fixedImage space
            regParams = calcImageTransform(sceneGeometryIn,x,cameraOffsetPoint);
            
            % Update the image            
            displayImage = updateFrame(displayImage,regParams,cameraOffsetPoint);
        end
            
        % Display the image
        imshow(displayImage,[],'Border','tight');
        ax = gca;
        ax.Toolbar = [];
        hold on
        % Add a label
        if showMoving
            text(20,30,'MOVING', 'Color', 'r','Fontsize',16);
        else
            text(20,30,'FIXED', 'Color', 'g','Fontsize',16);
        end
        
        % Add a marker for the camera CoP
        plot(cameraOffsetPoint(1),cameraOffsetPoint(2),'+c');
        
        keyAction = waitforbuttonpress;
        if keyAction
            keyChoiceValue = double(get(gcf,'CurrentCharacter'));
            switch keyChoiceValue
                case 28
                    x(2)=x(2)-0.1;
                case 29
                    x(2)=x(2)+0.1;
                case 30
                    x(3)=x(3)-0.1;
                case 31
                    x(3)=x(3)+0.1;
                case {45 95}
                    x(4)=x(4)-1;
                case {61 43}
                    x(4)=x(4)+1;
                case 102
                    showMoving = ~showMoving;
                case 106
                    x(1)=x(1)+0.5;
                case 107
                    x(1)=x(1)-0.5;
                case 112
                    showPerimeter = ~showPerimeter;
                case 109
                    showModel = ~showModel;
                case 27
                    text_str = 'finishing...';
                    addAnnotation(text_str);
                    stillWorking = false;
                otherwise
                    % Nothing to do
            end
            
        end
    end
    
    close(figHandle);
end


%% Create the adjusted sceneGeometry
sceneGeometryAdjusted = updateSceneGeometry( sceneGeometryIn, x );


%% Update the fixation angles
[~,modeledEyePose] = calcGlintGazeError( sceneGeometryAdjusted, args{:}, keyVals{:} );
medianEyePoseFixed = median(modeledEyePose);
sceneGeometryAdjusted.screenPosition.fixationEyePose = medianEyePoseFixed(1:2)';


%% Create and save a diagnostic figure
if p.Results.saveDiagnosticPlot
    
    % Moving frame -- Typically the gazeCal source
    displayImage = videoFrameMoving;
    tmpFig = figure('visible','off');
    renderEyePose(eyePoseMoving, sceneGeometryIn, ...
        'newFigure', false, 'visible', false, ...
        'backgroundImage',displayImage, ...
        'showAzimuthPlane', true, ...
                'modelEyeLabelNames', {'retina' 'pupilEllipse' 'cornea' 'glint_01'}, ...
                'modelEyePlotColors', {'.w' '-g' '.y' 'xr'}, ...
                'modelEyeSymbolSizeScaler',1.5,...
                'modelEyeAlpha', [0.25 0.25 0.25 1]);            
    text(20,30,sceneGeometryInStem, 'Color', 'g','Fontsize',16,'Interpreter','none');
    msg = ['frame ' num2str(referenceFrameMoving)];
    addAnnotation(msg);
    hold on
    plot([size(displayImage,2)/2, size(displayImage,2)/2],[0 size(displayImage,2)],'-b');
    plot([0 size(displayImage,1)],[size(displayImage,1)/2, size(displayImage,1)/2],'-b');
    tmpFrame = getframe(gcf);
    imageSet(1) = {tmpFrame.cdata};
    close(tmpFig);

        
    % Fixed frame -- The acquisition for which we have new sceneGeometry
    displayImage = videoFrameFixed;
    eyePoseFixed = eyePoseEllipseFit(XpFixed, YpFixed, sceneGeometryAdjusted,'glintCoord',glintCoordFixed,keyVals{:});
    tmpFig = figure('visible','off');
    renderEyePose(eyePoseFixed, sceneGeometryAdjusted, ...
        'newFigure', false, 'visible', false, ...
        'backgroundImage',displayImage, ...
        'showAzimuthPlane', true, ...
                'modelEyeLabelNames', {'retina' 'pupilEllipse' 'cornea' 'glint_01'}, ...
                'modelEyePlotColors', {'.w' '-g' '.y' 'xr'}, ...
                'modelEyeSymbolSizeScaler',1.5,...
                'modelEyeAlpha', [0.25 0.25 0.25 1]);            
    text(20,30,sceneGeometryOutStem, 'Color', 'r','Fontsize',16,'Interpreter','none');
    msg = ['frame ' num2str(referenceFrameFixed)];
    addAnnotation(msg);
    % Add cross hairs
    hold on
    plot([size(displayImage,2)/2, size(displayImage,2)/2],[0 size(displayImage,2)],'-b');
    plot([0 size(displayImage,1)],[size(displayImage,1)/2, size(displayImage,1)/2],'-b');
    tmpFrame = getframe(gcf);
    imageSet(2) = {tmpFrame.cdata};
    close(tmpFig);
    
    
    % Difference image
    regParams = calcImageTransform(sceneGeometryIn,x,cameraOffsetPoint);
    adjMovingFrame = updateFrame(videoFrameMoving,regParams,cameraOffsetPoint);
    displayImage = videoFrameFixed - double(adjMovingFrame);
    tmpFig = figure('visible','off');
    imshow(displayImage,[], 'Border', 'tight');
    text(20,30,'Difference', 'Color', 'w','Fontsize',16,'Interpreter','none');
    tmpFrame = getframe(gcf);
    imageSet(3) = {tmpFrame.cdata};
    close(tmpFig);
    
    % Prepare the figure
    figHandle=figure('visible','off');
    set(gcf,'PaperOrientation','portrait');
    
    set(figHandle, 'Units','inches')
    height = 4;
    width = 12;
    
    % The last two parameters of 'Position' define the figure size
    set(figHandle, 'Position',[25 5 width height],...
        'PaperSize',[width height],...
        'PaperPositionMode','auto',...
        'Color','w',...
        'Renderer','painters'...
        );
    
    % Post the montage of the imageSet
    montage(imageSet,'Size', [1 3]);
    
    % Post the title
    pathParts = strsplit(sceneGeometryInPath,filesep);
    titleString = [fullfile(pathParts{end-4:end-2})];
    title(titleString,'Interpreter','none')
    
    % Report the alignment method
    annotation('textbox', [0.15, .125, 0, 0], 'string', alignMethod,'FontWeight','bold','FitBoxToText','on','LineStyle','none','HorizontalAlignment','left','Interpreter','none')     
    
    % Report the run lengths
    msg = sprintf('runLength fixed = %2.0f',runLengthFixed);
    annotation('textbox', [0.75, .125, 0, 0], 'string', msg,'FitBoxToText','on','LineStyle','none','HorizontalAlignment','left','Interpreter','none')     
    
    % Add a text summary below. If any delta fixation angle is geater than
    % 1 deg, print the message text in red to alert that this was a large
    % eye rotation change.
    deltaX = x-x0;
    deltaPose = medianEyePoseFixed - eyePoseMoving;
    msg = sprintf('delta torsion [deg] = %2.1f',deltaX(1));
    annotation('textbox', [0.5, .175, 0, 0], 'string', msg,'FitBoxToText','on','LineStyle','none','HorizontalAlignment','center','Interpreter','none')
    msg = sprintf('delta translation [mm] [x; y; z] = [%2.3f; %2.3f; %2.3f]',deltaX(2:4));
    annotation('textbox', [0.5, .125, 0, 0], 'string', msg,'FitBoxToText','on','LineStyle','none','HorizontalAlignment','center','Interpreter','none')
    msg = sprintf('delta eye pose [azi, ele, tor, radius] = [%2.3f, %2.3f, %2.3f, %2.3f]',deltaPose);
    msgColor = 'black';
    if any(abs(deltaPose(1:3)) > 0) 
        msgColor = 'red';
    end
    annotation('textbox', [0.5, .075, 0, 0], 'string', msg,'Color',msgColor,'FitBoxToText','on','LineStyle','none','HorizontalAlignment','center','Interpreter','none')
    
    % Save and close the figure
    tmp = fullfile(sceneGeometryOutPath,[sceneGeometryOutStem '_sceneSync_QA.pdf']);
    print(figHandle,tmp,'-dpdf','-bestfit');
    close(figHandle);
    
end


%% Save the adjusted sceneGeometry
if p.Results.saveAdjustedSceneGeometry
    
    % If the movingFrame is from a point after the start of the
    % acquisition, adjust the camera translation position once more to
    % account for any head motion that took place between the start of the
    % acquisition and the reference frame
    sceneGeometryAdjusted.cameraPosition.translation = ...
        sceneGeometryAdjusted.cameraPosition.translation + ...
        relativeCameraPosition.values(:,referenceFrameFixed);
    
    % Add the meta data
    sceneGeometryAdjusted.meta.syncSceneGeometry = p.Results;
    sceneGeometryAdjusted.meta.syncSceneGeometry.x = x;
    
    % Set the variable name
    sceneGeometry = sceneGeometryAdjusted;
    
    % Save it
    tmp = fullfile(sceneGeometryOutPath,[sceneGeometryOutStem '_sceneGeometry.mat']);
    save(tmp,'sceneGeometry');
end

% If we are in display mode, report the values to the console
if p.Results.displayMode
    tmp=strsplit(sceneGeometryInPath,filesep);
    outline = [tmp{end-3} char(9) tmp{end-2} char(9) 'fixed: ' sceneGeometryInStem ', moving: ' sceneGeometryOutStem '\n'] ;
    fprintf(outline)
    outline = ['alignMethod' char(9) 'x' char(9) 'eyePose\n'];
    fprintf(outline)
    outline = sprintf(['{''' alignMethod '''}' char(9) '[ %2.2f, %2.2f, %2.2f, %2.2f, %2.2f, %2.2f, %2.2f, %2.2f ]' char(9) '[ %2.2f, %2.2f, %2.2f, %2.2f ]\n'],x,medianEyePoseFixed);
    fprintf(outline)
    fprintf('\n')
end


end % Main function





%% LOCAL FUNCTIONS


function regParams = calcImageTransform(sceneGeometry,x,cameraOffsetPoint)
% Determine the rotation and translation matrices that describe the change
% in an image induced by the updated sceneParameters
f = updateSceneGeometry(sceneGeometry,x);

pupilEllipseA1 = projectModelEye([ 0 1 0 1],sceneGeometry,'pupilRayFunc',[]);
pupilEllipseA2 = projectModelEye([-1 0 0 1],sceneGeometry,'pupilRayFunc',[]);
pupilEllipseA3 = projectModelEye([ 1 0 0 1],sceneGeometry,'pupilRayFunc',[]);
pupilEllipseB1 = projectModelEye([ 0 1 0 1],f,'pupilRayFunc',[]);
pupilEllipseB2 = projectModelEye([-1 0 0 1],f,'pupilRayFunc',[]);
pupilEllipseB3 = projectModelEye([ 1 0 0 1],f,'pupilRayFunc',[]);

A = [pupilEllipseA1(1:2)',pupilEllipseA2(1:2)',pupilEllipseA3(1:2)']-cameraOffsetPoint';
B = [pupilEllipseB1(1:2)',pupilEllipseB2(1:2)',pupilEllipseB3(1:2)']-cameraOffsetPoint';

regParams = absor(...
    A,...
    B,...
    'doScale',true,...
    'doTrans',true);
end


function displayImage = updateFrame(movingFrame,regParams,cameraOffsetPoint)

for ii=1:size(movingFrame,3)
    tmpImage = squeeze(movingFrame(:,:,ii));
    % Embed the movingFrame within a larger image that is padded
    % with mid-point background values
    padVals = round(size(tmpImage)./2);
    tmpImagePad = zeros(size(tmpImage)+padVals.*2)+125;
    tmpImagePad(padVals(1)+1:padVals(1)+size(tmpImage,1), ...
        padVals(2)+1:padVals(2)+size(tmpImage,2) ) = tmpImage;
    tmpImage = tmpImagePad;
    % Rotate the image
    tmpImage = imrotateAround(tmpImage, cameraOffsetPoint(2), cameraOffsetPoint(1), regParams.theta, 'bicubic');
    % Apply the x and y translation
    tmpImage = imtranslate(tmpImage,regParams.t','method','cubic');
    % Apply the scaling
    tmpImage = HardZoom(tmpImage,regParams.s);
    % Crop out the padding
    tmpImage = tmpImage(padVals(1)+1:padVals(1)+size(movingFrame,1), ...
        padVals(2)+1:padVals(2)+size(movingFrame,2));
    % Store this dimension
    displayImage(:,:,ii)=uint8(tmpImage);
end

end



function annotHandle = addAnnotation(text_str)
annotHandle = annotation('textbox',...
    [.80 .85 .1 .1],...
    'HorizontalAlignment','center',...
    'VerticalAlignment','middle',...
    'Margin',1,...
    'String',text_str,...
    'FontSize',9,...
    'FontName','Helvetica',...
    'EdgeColor',[1 1 1],...
    'LineWidth',1,...
    'BackgroundColor',[0.9  0.9 0.9],...
    'Color',[1 0 0]);
drawnow
end


function OutPicture = HardZoom(InPicture, ZoomFactor)
      % Si el factor de escala es 1, no se hace nada
      if ZoomFactor == 1
          OutPicture = InPicture;
          return;
      end
      % Se obtienen las dimensiones de las imágenes
      ySize = size(InPicture, 1);
      xSize = size(InPicture, 2);
      zSize = size(InPicture, 3);
      yCrop = floor(ySize / 2 * abs(ZoomFactor - 1));
      xCrop = floor(xSize / 2 * abs(ZoomFactor - 1));
      % Si el factor de escala es 0 se devuelve una imagen en negro
      if ZoomFactor == 0
          OutPicture = uint8(zeros(ySize, xSize, zSize));
          return;
      end
      % Se reescala y se reposiciona en en centro
      zoomPicture = imresize(InPicture, ZoomFactor);
      ySizeZ = size(zoomPicture, 1);
      xSizeZ = size(zoomPicture, 2);      
      if ZoomFactor > 1
          OutPicture = zoomPicture( 1+yCrop:yCrop+ySize, 1+xCrop:xCrop+xSize, :);
      else
          OutPicture = uint8(zeros(ySize, xSize, zSize));
          OutPicture( 1+yCrop:yCrop+ySizeZ, 1+xCrop:xCrop+xSizeZ, :) = zoomPicture;
      end
end


function output = imrotateAround(image, pointY, pointX, angle, varargin)
% ROTATEAROUND rotates an image.
%   ROTATED=ROTATEAROUND(IMAGE, POINTY, POINTX, ANGLE) rotates IMAGE around
%   the point [POINTY, POINTX] by ANGLE degrees. To rotate the image
%   clockwise, specify a negative value for ANGLE.
%
%   ROTATED=ROTATEAROUND(IMAGE, POINTY, POINTX, ANGLE, METHOD) rotates the
%   image with specified method:
%       'nearest'       Nearest-neighbor interpolation
%       'bilinear'      Bilinear interpolation
%       'bicubic'       Bicubic interpolation
%    The default is fast 'nearest'. Switch to 'bicubic' for nicer results.
%
%   Example
%   -------
%       imshow(rotateAround(imread('eight.tif'), 1, 1, 10));
%
%   See also IMROTATE, PADARRAY.
%   Contributed by Jan Motl (jan@motl.us)
%   $Revision: 1.2 $  $Date: 2014/05/01 12:08:01 $
% Parameter checking.
numvarargs = length(varargin);
if numvarargs > 1
    error('myfuns:somefun2Alt:TooManyInputs', ...
        'requires at most 1 optional input');
end
optargs = {'nearest'};    % Set defaults for optional inputs
optargs(1:numvarargs) = varargin;
[method] = optargs{:};    % Place optional args in memorable variable names
% Initialization.
[imageHeight, imageWidth, ~] = size(image);
centerX = floor(imageWidth/2+1);
centerY = floor(imageHeight/2+1);
dy = centerY-pointY;
dx = centerX-pointX;
% How much would the "rotate around" point shift if the
% image was rotated about the image center.
[theta, rho] = cart2pol(-dx,dy);
[newX, newY] = pol2cart(theta+angle*(pi/180), rho);
shiftX = round(pointX-(centerX+newX));
shiftY = round(pointY-(centerY-newY));
% Pad the image to preserve the whole image during the rotation.
padX = imageHeight;
padY = imageWidth;
padded = padarray(image, [padY padX],125);
% Rotate the image around the center.
rot = imrotate(padded, angle, method, 'crop');
% Crop the image.
output = rot(padY+1-shiftY:end-padY-shiftY, padX+1-shiftX:end-padX-shiftX, :);
end


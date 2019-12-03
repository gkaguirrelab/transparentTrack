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
%   sceneGeometryFileNameToSync - Full path to the .mat file that contains
%                           the sceneGeometry to be used.
%
% Examples:
%{
    % Invoke the file picker GUI
    syncSceneGeometry('','displayMode',true,'alignMethod','shape');
%}


%% Parse vargin for options passed here
p = inputParser; p.KeepUnmatched = true;

% Required
p.addRequired('pupilFileName',@ischar);

% Optional display and I/O params
p.addParameter('sceneGeometryFileNameToSync','',@(x)(ischar(x) | iscell(x)));
p.addParameter('verbose',true,@islogical);
p.addParameter('displayMode',false,@islogical);
p.addParameter('saveAdjustedSceneGeometry',true,@islogical);
p.addParameter('saveDiagnosticPlot',true,@islogical);
p.addParameter('doNotSyncSceneToItself',true,@islogical);

% Optional fitting params
p.addParameter('alignMethod','shape',@(x)(ischar(x) | iscell(x)));
p.addParameter('deltaPix',[],@isnumeric);
p.addParameter('deltaDeg',[],@isnumeric);
p.addParameter('deltaScale',[],@isnumeric);
p.addParameter('deltaPose',[],@isnumeric);
p.addParameter('eyePositionTargetLengthFrames',30,@isscalar);
p.addParameter('gazeErrorThreshTol',0.25,@isscalar);


%% Parse and check the parameters
p.parse(pupilFileName, varargin{:});


%% Load sceneGeometryIn files
% This sceneGeometry--and its associated video, perimeter, and pupil
% data--constitute the "fixed" measurements.

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
% video frame measurements tht best represent the appearance of the eye at
% this position.

% Get the camera offset point. This is where the center of the video image
% lands on the camera sensor array.
cameraOffsetPoint = [sceneGeometryIn.cameraIntrinsic.matrix(1,3), ...
    sceneGeometryIn.cameraIntrinsic.matrix(2,3)];

% Obtain the eye rotation values from the pupilData, and convert these into
% gaze position on the screen.
tmp = pupilData.radiusSmoothed.eyePoses.values(:,1:2)';
gazePosition = (sceneGeometryIn.screenPosition.R * tmp + sceneGeometryIn.screenPosition.fixationAngles(1:2)')';

% As we will be trying to find the frames during which the gaze was
% directed towards the center of the screen, non-zero gaze position counts
% as error in the search below.
eyeMatchError = sqrt(sum(gazePosition.^2,2));


%% Find the set of frames for the fixed image

% Anonymous function to grab the indicies of when runs of frames begin
pullStartIndices = @(vec) vec(1:2:end-1);

% Anonymous function to grab the length of each run
pullRunLengths = @(vec) vec(2:2:end)-pullStartIndices(vec);

% Find the minimum fixation error threshold that results in a run of
% consecutive frames at fixation of the target length.
targetLength = p.Results.eyePositionTargetLengthFrames;

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
    
    % The minimum threshold eyeMatchError that results in a run length that
    % matches the target run length.
    threshValFixed = fzero(myObj,0.5,options);
    
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
    warning('Unable to find a suitable set of frames from sceneGeometryIn')
    return
end

% Find the start point of this run of frames
runLengths = pullRunLengths(runStarts(threshValFixed));
runIndices = pullStartIndices(runStarts(threshValFixed));
runLengthFixed = targetLength-myObj(threshValFixed);
startIndexFixed = runIndices(runLengths == runLengthFixed);
startIndexFixed = startIndexFixed(1);

% Check that the signed gaze error is acceptable. If not, adopt an
% alternate strategy for picking the frames
signedError = sqrt(sum(gazePosition(startIndexFixed:startIndexFixed+runLengthFixed-1,1))^2 + ...
    sum(gazePosition(startIndexFixed:startIndexFixed+runLengthFixed-1,2))^2);
if signedError > p.Results.gazeErrorThreshTol
    % Functions to return signed and unsigned error as a function of
    % runLength
    signedErrorByRunLength = @(k) sqrt(movsum(gazePosition(:,1),k,'Endpoints','fill').^2 + movsum(gazePosition(:,2),k,'Endpoints','fill').^2);
    unsignedErrorByRunLength = @(k) movsum( sqrt(gazePosition(:,1).^2 + gazePosition(:,2).^2),k,'Endpoints','fill' )./k;
    % Reset the targetLength
    targetLength = p.Results.eyePositionTargetLengthFrames;
    % Loop through and find the maximum target length that produces an
    % acceptable signed and unsigned error
    stillWorking = true;
    while stillWorking
        validIdx = find ( double( signedErrorByRunLength(targetLength) < p.Results.gazeErrorThreshTol ) .* ...
            double( unsignedErrorByRunLength(targetLength) < p.Results.gazeErrorThreshTol ) );        
        if length(validIdx) > 1
            stillWorking = false;
            startIndexFixed = validIdx(1) - floor(targetLength/2);
            runLengthFixed = targetLength;
        else
            targetLength = targetLength - 1;
        end
        % if targetLength has hit 1, we have run out of chances. Just take
        % the frame with the least signed error
        if targetLength == 1
            [~,startIndexFixed] = min(signedErrorByRunLength(targetLength));
            runLengthFixed = 1;
            stillWorking = false;
        end
    end
end

% Find the frame with the lowest ellipse RMSE during this period. We will
% desginate this frame the referenceFrame for the fixed (sceneGeometryIn)
% measurements
rmseVals = pupilData.radiusSmoothed.ellipses.RMSE(startIndexFixed:startIndexFixed+runLengthFixed-1);
referenceFrameFixed = startIndexFixed + find(rmseVals == min(rmseVals)) - 1;

% Obtain and store the pupil perimeter points for this frame
XpFixed = perimeter.data{referenceFrameFixed}.Xp;
YpFixed = perimeter.data{referenceFrameFixed}.Yp;

% The eyePose for this reference frame is set to the fixation value for the
% sceneGeometry, with the pupil size set to the median value for this run
% of frames.
eyePoseFixed = nanmedian(pupilData.radiusSmoothed.eyePoses.values(startIndexFixed:startIndexFixed+runLengthFixed-1,:));
eyePoseFixed(1:3) = -sceneGeometryIn.screenPosition.fixationAngles;

% Load in the median image from the period of fixation for sceneGeometryIn.
% This is the "fixed" frame.
tmp = fullfile(sceneGeometryInPath,[sceneGeometryInStem '_gray.avi']);
videoFrameFixed = makeMedianVideoImage(tmp,'startFrame',startIndexFixed,'nFrames',runLengthFixed);

% Find the median pupil center across the frames for the run of frames
pupilCenterPixelsFixed = [ ...
    nanmedian(pupilData.radiusSmoothed.ellipses.values(startIndexFixed:startIndexFixed+runLengthFixed-1,1)), ...
    nanmedian(pupilData.radiusSmoothed.ellipses.values(startIndexFixed:startIndexFixed+runLengthFixed-1,2)) ];

% Find the median shape of the pupil during the run of frames, expressed as
% theta and rho values (SEE: csaEllipseError)
pupilRhoShapeFixed = nanmedian(pupilData.radiusSmoothed.ellipses.values(startIndexFixed:startIndexFixed+runLengthFixed-1,4));
pupilRhoShapeFixed = 1-sqrt(1-pupilRhoShapeFixed^2);
pupilThetaShapeFixed = nanmedian(pupilData.radiusSmoothed.ellipses.values(startIndexFixed:startIndexFixed+runLengthFixed-1,5));
pupilThetaShapeFixed = pupilThetaShapeFixed*2;


%% Load sceneGeometryOut files
% This acquisition--and its associated video, perimeter, and pupil
% data--constitute the "moving" measurements. Our goal is to adjust the
% sceneGeometryIn to best match this acquisition.

% If the pupilFileName is not defined, offer some choices
if isempty(pupilFileName)
    % Get a list of all gray.avi videos in this directory
    fileList = dir(fullfile(sceneGeometryInPath,'*_pupil.mat'));
    
    % Exclude the video that is the source of the fixed image
    keep=cellfun(@(x) ~strcmp(x,[sceneGeometryInStem '_pupil.mat']),extractfield(fileList,'name'));
    fileList = fileList(keep);
    
    % Ask the operator which of the videos we wish to adjust
    fprintf('\n\nSelect the pupil data to adjust:\n')
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
tmp =  fullfile(sceneGeometryOutPath,[sceneGeometryOutStem '_correctedPerimeter.mat']);
load(tmp,'perimeter');
tmp =  fullfile(sceneGeometryOutPath,[sceneGeometryOutStem '_relativeCameraPosition.mat']);
load(tmp,'relativeCameraPosition');

% Identify the acqStartTimeMoving, which is the time point at which the
% fMRI acquisition began
[~, acqStartFrameMoving] = min(abs(timebase.values));


%% Check if we are syncing a sceneGeometry to itself
if p.Results.doNotSyncSceneToItself && strcmp(sceneGeometryInStem,sceneGeometryOutStem)
    if p.Results.verbose
        fprintf('Detected that the sceneGeometry source is the same as the acquisition; returning.\n');
    end
    return
end


%% Find the target window for the acqusition
% Identify a window of frames from the acquisition during which the eye is
% fixating upon "screen position" [0, 0], and an error metric for the
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
        windowEnd = acqStartFrameMoving;
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
        % a movie is the center of the screen.
        % Find the period after to the start of the scan when the eye was
        % in the most consistent position, and closest to the median
        % position.
        windowStart = acqStartFrameMoving;
        windowEnd = acqStartFrameMoving+600;
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
        windowStart = acqStartFrameMoving;
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

% Anonynous function that provides the lengths of runs of frames for which
% the eyeMatch error is below a threshold.
runStarts = @(thresh) find(diff([0,(eyeMatchError < thresh)',0]==1));

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
    threshValMoving = fzero(myObj,x0,options);
    
    % Check to see if the search has met our threshold criteria, or if we
    % have run the targetLength down as short as it can go.
    if threshValMoving < p.Results.gazeErrorThreshTol || targetLength == 1
        stillWorking = false;
    else
        targetLength = targetLength-1;
    end
end

% Check if we found a solution
if ~isfinite(threshValMoving)
    warning('Unable to find a suitable set of frames from the acquisition')
    return
end

% Find the start point of this run of frames
runLengths = pullRunLengths(runStarts(threshValMoving));
runIndices = pullStartIndices(runStarts(threshValMoving))+windowStart-1;
runLengthMoving = targetLength-myObj(threshValMoving);
startIndexMoving = runIndices(runLengths == runLengthMoving);
startIndexMoving = startIndexMoving(1);

% Find the frame with the lowest ellipse RMSE during the target window
rmseVals = pupilData.initial.ellipses.RMSE(startIndexMoving:startIndexMoving+runLengthMoving-1);
referenceFrameMoving = startIndexMoving + find(rmseVals == min(rmseVals)) - 1;


%% Derive properties from the acquisition
% Obtain some measurements of things from within the target window of the
% acqusition

% The pupil perimeter for the reference frame
XpMoving = perimeter.data{referenceFrameMoving}.Xp;
YpMoving = perimeter.data{referenceFrameMoving}.Yp;

% Obtain the median [x y] position of the pupil center during the target
% period of the moving image, and use this to determine the displacement
% (in pixels) from the [x y] position of the pupil center during the
% corresponding target period from the fixed (sceneGeometryIn) image.

% Get the pupil center for the fames from the moving video
pupilCenterPixelsMoving = [ ...
    nanmedian(pupilData.initial.ellipses.values(startIndexMoving:startIndexMoving+runLengthMoving-1,1)), ...
    nanmedian(pupilData.initial.ellipses.values(startIndexMoving:startIndexMoving+runLengthMoving-1,2)) ];

% Load in the median image from the target period for acquisition.
% This is the "moving" frame.
tmp = fullfile(sceneGeometryOutPath,[sceneGeometryOutStem '_gray.avi']);
movingFrame = makeMedianVideoImage(tmp,'startFrame',startIndexMoving,'nFrames',runLengthMoving);


%% Set the displacement of the moving image in pixels and degrees
% This can either be passed as a key-value, or derived from the
% sceneGeometry and acquisition

if isempty(p.Results.deltaPix)
    % The adjustment is the difference in pupil centers from the fixed
    % and moving videos
    deltaPix = pupilCenterPixelsFixed - pupilCenterPixelsMoving;
else
    deltaPix = p.Results.deltaPix;
end

if isempty(p.Results.deltaScale)
    deltaScale = 1;
end


if isempty(p.Results.deltaDeg)
    % No change is made to the torsion unless we are in display mode
    deltaDeg = 0;
else
    deltaDeg = p.Results.deltaDeg;
end


%% Optional display mode
% If we are in display mode, offer an interface for the user to manually
% adjust the delta variables

if p.Results.displayMode
    
    % Provide some instructions for the operator
    fprintf('Adjust horizontal and vertical camera translation with the arrow keys.\n');
    fprintf('Adjust depth camera translation with + and -.\n');
    fprintf('Adjust camera torsion with j and k.\n');
    fprintf('Switch between moving and fixed image by pressing a.\n');
    fprintf('Turn on and off perimeter display with p.\n');
    fprintf('Turn on and off model display with m.\n');
    fprintf('Press esc to exit.\n\n');
    fprintf([sceneGeometryInPath '\n']);
    
    % Create a figure
    figHandle = figure();
    imshow(videoFrameFixed,[],'Border','tight');
    ax = gca;
    ax.Toolbar = [];
    hold on
    text(20,30,'FIXED', 'Color', 'g','Fontsize',16);
    
    % Prepare for the loop
    showMoving = true;
    showPerimeter=false;
    showModel=false;
    stillWorking = true;
    
    % Enter the while stillWorking loop
    while stillWorking
        
        % Prepare to update the image
        hold off
        
        if showMoving
            % Work with the moving frame
            displayImage = updateMovingFrame(movingFrame,deltaPix,deltaDeg,deltaScale,cameraOffsetPoint);
            
            % Update the perimeter points
            [XpDisplay, YpDisplay] = updatePerimeter(XpMoving,YpMoving,deltaPix,deltaDeg,deltaScale,cameraOffsetPoint);
            
            % Display the perimeter points
            if showPerimeter
                idx = sub2ind(size(displayImage),round(YpDisplay),round(XpDisplay));
                displayImage(idx)=255;
            end
            
            % Display the image
            imshow(displayImage,[],'Border','tight');
            ax = gca;
            ax.Toolbar = [];
            hold on
            text(20,30,'MOVING', 'Color', 'r','Fontsize',16);
            
        else
            % Work with the fixed frame
            displayImage = videoFrameFixed;
            
            % Update the perimeter
            XpDisplay = XpFixed; YpDisplay = YpFixed;
            
            % Display the perimeter points
            if showPerimeter
                idx = sub2ind(size(displayImage),round(YpDisplay),round(XpDisplay));
                displayImage(idx)=255;
            end
            
            imshow(displayImage,[],'Border','tight');
            ax = gca;
            ax.Toolbar = [];
            hold on
            text(20,30,'FIXED', 'Color', 'g','Fontsize',16);
        end
        
        % Show the eye model
        if showModel
            % Let the user know this will take a few seconds
            text_str = 'Updating model...';
            annotHandle = addAnnotation(text_str);
            % Obtain the eye pose from the adjusted perimeter
            if isempty(p.Results.deltaPose)
                eyePoseDisplay = eyePoseEllipseFit(XpDisplay, YpDisplay, ...
                    sceneGeometryIn,'x0',eyePoseFixed);
            else
                eyePoseDisplay = eyePoseFixed+p.Results.deltaPose;
            end
            % Render the eye model
            renderEyePose(eyePoseDisplay, sceneGeometryIn, ...
                'newFigure', false, 'visible', true, ...
                'showAzimuthPlane', true, ...
                'modelEyeLabelNames', {'retina' 'pupilEllipse' 'cornea'}, ...
                'modelEyePlotColors', {'.w' '-g' '.y'}, ...
                'modelEyeSymbolSizeScaler',1.5,...
                'modelEyeAlpha', 0.25);
            hold on
            % Remove the updating annotation
            delete(annotHandle);
        end
        
        % Add a marker for the camera CoP
        plot(cameraOffsetPoint(1),cameraOffsetPoint(2),'+c');
        
        keyAction = waitforbuttonpress;
        if keyAction
            keyChoiceValue = double(get(gcf,'CurrentCharacter'));
            switch keyChoiceValue
                case 28
                    deltaPix(1)=deltaPix(1)-1;
                case 29
                    deltaPix(1)=deltaPix(1)+1;
                case 30
                    deltaPix(2)=deltaPix(2)-1;
                case 31
                    deltaPix(2)=deltaPix(2)+1;
                case {45 95}
                    deltaScale = deltaScale+0.01;
                case {61 43}
                    deltaScale = deltaScale-0.01;
                case 97
                    showMoving = ~showMoving;
                case 106
                    deltaDeg = deltaDeg - 1;
                case 107
                    deltaDeg = deltaDeg + 1;
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
sceneGeometryAdjusted = sceneGeometryIn;

% Update the sceneGeometry torsion
sceneGeometryAdjusted.cameraPosition.torsion = sceneGeometryIn.cameraPosition.torsion - deltaDeg;

% Find the change in mm of extrinsic camera translation needed to shift the
% eye model the observed number of pixels
adjustedTranslation = calcCameraTranslationPixels(sceneGeometryAdjusted,eyePoseFixed,deltaPix);
deltaMM = sceneGeometryIn.cameraPosition.translation - adjustedTranslation;

% Update the sceneGeometry translation
sceneGeometryAdjusted.cameraPosition.translation = adjustedTranslation;

% Update the eye pose for the adjusted sceneGeometry
if ~isempty(p.Results.deltaPose)
	eyePoseAdjusted = eyePoseFixed+p.Results.deltaPose;
else
    for ii = 1:runLengthMoving
        % The pupil perimeter for the reference frame
        Xpt = perimeter.data{startIndexMoving+ii-1}.Xp;
        Ypt = perimeter.data{startIndexMoving+ii-1}.Yp;
        if ~isempty(Xpt)
            if length(Xpt)>5
                eyePoseByFrame(ii,:) = eyePoseEllipseFit(Xpt, Ypt, ...
                    sceneGeometryAdjusted,'x0',eyePoseFixed);
            else
                eyePoseByFrame(ii,:)=[nan nan nan nan];
            end
        else
            eyePoseByFrame(ii,:)=[nan nan nan nan];
        end
    end
    % Take the weighted median across eyePoses
    weights = weightFunc(pupilData,perimeter,startIndexMoving,runLengthMoving);
    nonNanFrames = ~isnan(eyePoseByFrame(:,1));
    for ii=1:4
        eyePoseAdjusted(ii) = medianw(eyePoseByFrame(nonNanFrames,ii),weights(nonNanFrames));
    end
end
sceneGeometryAdjusted.screenPosition.fixationAngles = -eyePoseAdjusted(1:3);
            

%% Create and save a diagnostic figure
if p.Results.saveDiagnosticPlot
    
    % Fixed frame
    displayImage = videoFrameFixed;
    tmpFig = figure('visible','off');
    renderEyePose(eyePoseFixed, sceneGeometryIn, ...
        'newFigure', false, 'visible', false, ...
        'backgroundImage',displayImage, ...
        'showAzimuthPlane', true, ...
        'modelEyeLabelNames', {'retina' 'pupilEllipse' 'cornea'}, ...
        'modelEyePlotColors', {'.w' '-g' '.y'}, ...
        'modelEyeSymbolSizeScaler',1.5,...
        'modelEyeAlpha', 0.25);
    text(20,30,sceneGeometryInStem, 'Color', 'r','Fontsize',16,'Interpreter','none');
    msg = ['frame ' num2str(referenceFrameFixed)];
    addAnnotation(msg);
    % Add cross hairs
    hold on
    plot([size(displayImage,2)/2, size(displayImage,2)/2],[0 size(displayImage,2)],'-b');
    plot([0 size(displayImage,1)],[size(displayImage,1)/2, size(displayImage,1)/2],'-b');
    tmpFrame = getframe(gcf);
    imageSet(1) = {tmpFrame.cdata};
    close(tmpFig);
    
    % Moving frame
    displayImage = movingFrame;
    tmpFig = figure('visible','off');
    renderEyePose(eyePoseAdjusted, sceneGeometryAdjusted, ...
        'newFigure', false, 'visible', false, ...
        'backgroundImage',displayImage, ...
        'showAzimuthPlane', true, ...
        'modelEyeLabelNames', {'retina' 'pupilEllipse' 'cornea'}, ...
        'modelEyePlotColors', {'.w' '-g' '.y'}, ...
        'modelEyeSymbolSizeScaler',1.5,...
        'modelEyeAlpha', 0.25);
    text(20,30,sceneGeometryOutStem, 'Color', 'g','Fontsize',16,'Interpreter','none');
    msg = ['frame ' num2str(referenceFrameMoving)];
    addAnnotation(msg);
    hold on
    plot([size(displayImage,2)/2, size(displayImage,2)/2],[0 size(displayImage,2)],'-b');
    plot([0 size(displayImage,1)],[size(displayImage,1)/2, size(displayImage,1)/2],'-b');
    tmpFrame = getframe(gcf);
    imageSet(2) = {tmpFrame.cdata};
    close(tmpFig);
    
    % Difference image
    adjMovingFrame = updateMovingFrame(movingFrame,deltaPix,deltaDeg,deltaScale,cameraOffsetPoint);
    displayImage = videoFrameFixed - adjMovingFrame;
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
    msg = sprintf('runLength [fixed, moving] = %2.0f, %2.0f',runLengthFixed,runLengthMoving);
    annotation('textbox', [0.75, .125, 0, 0], 'string', msg,'FitBoxToText','on','LineStyle','none','HorizontalAlignment','left','Interpreter','none')     
    
    % Add a text summary below. If any delta fixation angle is geater than
    % 1 deg, print the message text in red to alert that this was a large
    % eye rotation change.
    msg = sprintf('delta pixels = [x, y] = [%2.1f; %2.1f]',deltaPix);
    annotation('textbox', [0.5, .175, 0, 0], 'string', msg,'FitBoxToText','on','LineStyle','none','HorizontalAlignment','center','Interpreter','none')
    msg = sprintf('delta translation [mm] [x; y; z] = [%2.3f; %2.3f; %2.3f]',deltaMM);
    annotation('textbox', [0.5, .125, 0, 0], 'string', msg,'FitBoxToText','on','LineStyle','none','HorizontalAlignment','center','Interpreter','none')
    msg = sprintf('delta eye pose [azi, ele, tor, radius] = [%2.3f, %2.3f, %2.3f, %2.3f]',eyePoseAdjusted-eyePoseFixed);
    msgColor = 'black';
    if any(abs(eyePoseFixed(1:3)-eyePoseAdjusted(1:3)) > 0.5) 
        msgColor = '#ffa500'; % orange
    end
    if any(abs(eyePoseFixed(1:3)-eyePoseAdjusted(1:3)) > 1)
        msgColor = 'red';
    end
    annotation('textbox', [0.5, .075, 0, 0], 'string', msg,'Color',msgColor,'FitBoxToText','on','LineStyle','none','HorizontalAlignment','center','Interpreter','none')
    msg = sprintf('delta torsion [deg] = %2.3f',deltaDeg);
    annotation('textbox', [0.5, .025, 0, 0], 'string', msg,'FitBoxToText','on','LineStyle','none','HorizontalAlignment','center','Interpreter','none')
    
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
        relativeCameraPosition.values(:,referenceFrameMoving);
    
    % Add the meta data
    sceneGeometryAdjusted.meta.syncSceneGeometry = p;
    
    % Set the variable name
    sceneGeometry = sceneGeometryAdjusted;
    
    % Save it
    tmp = fullfile(sceneGeometryOutPath,[sceneGeometryOutStem '_sceneGeometry.mat']);
    save(tmp,'sceneGeometry');
end


end % Main function





%% LOCAL FUNCTIONS

function [Xp, Yp] = updatePerimeter(Xp,Yp,deltaPix,deltaDeg,deltaScale,cameraOffsetPoint)

% Create a matrix of the perimeter points
v = [Xp';Yp'];

% Create the translation matrix
t = repmat([deltaPix(1); deltaPix(2)], 1, length(Xp));

% Translate the points
v = v+t;

% Set up the rotation matrix
center = repmat([cameraOffsetPoint(1); cameraOffsetPoint(2)], 1, length(Xp));
theta = deg2rad(deltaDeg);
R = [cos(theta) -sin(theta); sin(theta) cos(theta)];

% Apply the rotation
v = R*(v - center) + center;

% Extract the Xp and Yp vectors
Xp = v(1,:)';
Yp = v(2,:)';

end


function p = calcCameraTranslationPixels(sceneGeometry,eyePose,deltaPix)
% Find the change in the extrinsic camera translation needed to shift
% the eye model the observed number of pixels for an eye with zero rotation
p0 = sceneGeometry.cameraPosition.translation;
ub = sceneGeometry.cameraPosition.translation + [20; 20; 0];
lb = sceneGeometry.cameraPosition.translation - [20; 20; 0];
place = {'cameraPosition' 'translation'};
mySG = @(p) setfield(sceneGeometry,place{:},p);
pupilCenter = @(k) k(1:2);
targetPupilCenter = pupilCenter(pupilProjection_fwd(eyePose,sceneGeometry)) - deltaPix;
myError = @(p) norm(targetPupilCenter-pupilCenter(pupilProjection_fwd(eyePose,mySG(p))));
options = optimoptions(@fmincon,'Diagnostics','off','Display','off');
p = fmincon(myError,p0,[],[],[],[],lb,ub,[],options);
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


function displayImage = updateMovingFrame(movingFrame,deltaPix,deltaDeg,deltaScale,cameraOffsetPoint)

% Embed the movingFrame within a larger image that is padded
% with mid-point background values
padVals = round(size(movingFrame)./2);
displayImagePad = zeros(size(movingFrame)+padVals.*2)+125;
displayImagePad(padVals(1)+1:padVals(1)+size(movingFrame,1), ...
    padVals(2)+1:padVals(2)+size(movingFrame,2) ) = movingFrame;
displayImage = displayImagePad;
% Apply the scaling
scaledImage = imresize(displayImage,deltaScale);
% Apply the x and y translation
displayImage = imtranslate(displayImage,deltaPix,'method','cubic');
% Crop out the padding
displayImage = displayImage(padVals(1)+1:padVals(1)+size(movingFrame,1), ...
    padVals(2)+1:padVals(2)+size(movingFrame,2));
% Rotate the image
displayImage = imrotateAround(displayImage, cameraOffsetPoint(2), cameraOffsetPoint(1), -deltaDeg, 'bicubic');

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


function weights = weightFunc(pupilData,perimeter,startIndex,runLength)
%% Calculate likelhood SD across frames
% Obtain a measure for each frame of how completely the perimeter points
% define a full, 360 degrees around the pupil. This index of coverage of
% the pupil perimeter is distVals. If there is a perfectly uniform angular
% distribution of points in space around the pupil perimeter, then the
% distVals value will be zero. If there are perimeter points only at a
% single angular location around the pupil cirle, then the distVal will 1.

% The likelihood SD is based upon the RMSE of the fit of the elipse to the
% perimeter points for each frame
RMSE = pupilData.initial.ellipses.RMSE(startIndex:startIndex+runLength-1)';

% Define the bins over which the distribution of perimeter angles will be
% evaluated. 20 bins works pretty well.
nDivisions = 20;
histBins = linspace(-pi,pi,nDivisions);

% Anonymous function returns the linear non-uniformity of a set of values,
% ranging from 0 when perfectly uniform to 1 when completely non-uniform.
nonUniformity = @(x) (sum(abs(x/sum(x)-mean(x/sum(x))))/2)/(1-1/length(x));

% Loop over frames. Frames which have no perimeter points will be given a
% distVal of NaN.
for ii = 1:runLength
    
    % Obtain the center of this fitted ellipse
    centerX = pupilData.initial.ellipses.values(startIndex+ii-1,1);
    centerY = pupilData.initial.ellipses.values(startIndex+ii-1,2);
    
    % Obtain the set of perimeter points
    Xp = perimeter.data{startIndex+ii-1}.Xp;
    Yp = perimeter.data{startIndex+ii-1}.Xp;
    
    % Calculate the deviation of the distribution of points from uniform
    linearNonUniformity(ii) = nonUniformity(histcounts(atan2(Yp-centerY,Xp-centerX),histBins));
end

% Subject the linearNonUniformity vector to a non-linear transformation.
% This has the effect of changing 0 --> 0.1, 0.8 --> 1, and values > 0.8
% --> infinity. This causes pupil perimeters with support at a single
% location (as opposed to fully around the perimeter) to have a markedly
% increased likelihood SD. Also, set InF values to something arbitrarily
% large.
distVals = (1./(1-sqrt(linearNonUniformity)))./10;
distVals(isinf(distVals)) = 1e20;

% The weight for each frame is the inverse of the RMSE multiplied by the distVal
weights = 1./(distVals.*RMSE);
end

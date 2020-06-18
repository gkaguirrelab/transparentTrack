function saveEyeModelMontage(obj,fileNameSuffix)
% Saves a montage of video frames and the superimposed model eye
%
% Syntax:
%  obj.saveEyeModelMontage(fileNameSuffix)
%
% Description:
%   Saves a PNG image of the posed model eye superimposed on the video
%   frames that were the source of the scene information.
%
% Inputs:
%   fileNameSuffix        - Char vector. A suffix for the saved file name.
%                           The calling function might pass '1', '2', etc
%                           to label different stages of the evolution of
%                           the scene object.
%
% Outputs:
%   none
%


%% Setup variables

% Obtain variables from the object
sceneGeometry = obj.sceneGeometry;
modelEyePose = obj.modelEyePose;
modelCameraTrans = obj.modelCameraTrans;
frameSet = obj.frameSet;
videoStemName = obj.videoStemName;
perimeter = obj.perimeter;
glintDataX = obj.glintDataX;
glintDataY = obj.glintDataY;
gazeTargets = obj.gazeTargets;

% The name of the video file source, and the name of the montage figure to
% be saved.
grayVideoName = [videoStemName '_gray.avi'];
montageFileName = [videoStemName '_sceneGeometry_eyeModelMontage' fileNameSuffix '.png'];


%% Exit early if no video or frames to show
if ~exist(grayVideoName,'file') || isempty(frameSet)
    return
end


%% Calculate montage positions
% We will create a square grid of frames, and try to position the frames in
% the grid by their eye pose.

% The number of grid squares
nFrames = length(frameSet);
nGrid = ceil(sqrt(nFrames));

% Convert the gazeTargets array into an ordered list to store the frames
pos = gazeTargets ./ max(unique(abs(gazeTargets)));

% Reverse the sign of the elevation row, as positive eye movements
% correspond to the eye moving upwards in the image
pos(2,:) = -pos(2,:);

% Convet the positions into an ordered list
montageOrder = (pos(2,:)+nGrid-1).*nGrid+pos(1,:)-1;

% Sort the ellipse array list so that the frames appear in temporal order
[frameSet, sortOrder] = sort(frameSet);
modelEyePose = modelEyePose(sortOrder,:);
modelCameraTrans = modelCameraTrans(:,sortOrder);
perimeter = perimeter(sortOrder);
glintDataX = glintDataX(sortOrder);
glintDataY = glintDataY(sortOrder);

% If we have a full supply of gazeTargets, use them to define the montage
% order. Otherwise, just go with the sortOrder
if ~any(isnan(sum(gazeTargets))) && ~any(isnan(sum(pos)))
    montageOrder = montageOrder(sortOrder);
else
    montageOrder = sortOrder;
end


%% Prepare the video object and figure

% Open the video object
videoInObj = videoIOWrapper(grayVideoName,'ioAction','read');

% Get the video properties
videoSizeX = videoInObj.Width;
videoSizeY = videoInObj.Height;

% Define a variable to hold the selected frames
framesToMontage = zeros(videoSizeY,videoSizeX,3,nGrid^2,'uint8');

% Define a non-visible figure
figHandle = figure( 'Visible', 'off');
hAxes = gca();

% A counter for which frame of the video is currently available to read
curFrame = 0;


%% Display the frames
% We will read through every frame of the video, and pause for plotting
% when we hit a frame on our list to be displayed. This inefficient
% behavior is mandated by MATLAB not providing a method for random access
% of a video by frame, as opposed to an imprecise access by time stamp.

% Silence some errors that can arise during the rendering of the model eye,
% and during the montaging of the set of images
warningState = warning;
warning('off','projectModelEye:ellipseFitFailed');
warning('off','images:initSize:adjustingMag');

% Loop through the frameSet
for ii = 1:length(frameSet)
    
    % This is the next frame we would like to display
    idx = frameSet(ii);
    
    % Read up to the desired frame in the video
    while curFrame < idx
        sourceFrame = readFrame(videoInObj);
        curFrame = curFrame + 1;
    end
    
    % Display this frame and clean up the image
    imshow(sourceFrame,'Border', 'tight','Parent',hAxes);
    hold on
    axis off;
    
    % Add the pupil perimeter points
    plot(perimeter{ii}.Xp ,perimeter{ii}.Yp, '.w', 'MarkerSize', 1);
    
    % Add the rendered eye model
    eyePose = modelEyePose(ii,:);
    if ~any(isnan(eyePose))
        renderEyePose(eyePose, sceneGeometry, 'newFigure', false, ...
            'cameraTrans', modelCameraTrans(:,ii), ...
            'modelEyeLabelNames', {'retina' 'irisPerimeter' 'pupilEllipse' 'cornea' 'glint_01' 'glint_02' 'medialCanthus' 'lateralCanthus'}, ...
            'modelEyePlotColors', {'.w' '.b' '-g' '.y' 'xr' 'xr' 'Ok' 'Ok'}, ...
            'modelEyeAlpha', [0.25 0.25 0.25 0.25 1 1],...
            'modelEyeSymbolSizeScaler',1.5,...
            'showAzimuthPlane',true);
    end
    
    % Add the observed glint
    hold on
    plot(glintDataX(ii), glintDataY(ii),'or');
    
    % Get the frame
    drawnow;
    thisFrame=getframe(figHandle);
    
    % Add a text label for the frame number
    frameLabel = sprintf('frame: %d',idx);
    thisFrame.cdata = insertText(thisFrame.cdata,[20 20],frameLabel,'FontSize',30);
    
    % Store the frame if it is a valid size. It might be invalid if there
    % was something wrong with the rendering, for example.
    if all(size(squeeze(framesToMontage(:,:,:,ii)))==size(thisFrame.cdata))
        framesToMontage(:,:,:,montageOrder(ii)) = thisFrame.cdata;
    end
    
    % hold off
    hold off
end

% Close the temporary figure we used to render each frame
close(figHandle);

% Prepare a new, non-visible figure
figHandle = figure('Visible','off');
set(gcf,'PaperOrientation','landscape');
set(figHandle, 'Units','inches')
height = 6;
width = 11;

% The last two parameters of 'Position' define the figure size
set(figHandle, 'Position',[25 5 width height],...
    'PaperSize',[width height],...
    'PaperPositionMode','auto',...
    'Color','w');

% Create the montage
montage(framesToMontage);

% Save the montage
saveas(figHandle,montageFileName)

% Close the figure
close(figHandle)

% If Matlab is older than 2020a, we have to rotate the figure by 90 degrees
% clockwise, because I can't otherwise get the MATLAB plotting routines to
% output the image how I want it.
if verLessThan('matlab','9.8')
    A = imread(montageFileName);
    A = rot90(A,3);
    imwrite(A,montageFileName);
end

% close the video object
clear videoInObj

% Restore the warning state
warning(warningState);


end % saveEyeModelMontage

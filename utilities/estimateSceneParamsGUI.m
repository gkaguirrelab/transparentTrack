function [ initialParams ] = estimateSceneParamsGUI(sceneGeometryFileName)
% Adjust scene parameter values
%
% Syntax:
%  [ initialParams ] = estimatePipelineParamsGUI(grayVideoName, protocol)
%
% Description:



%% Input parser
p = inputParser; p.KeepUnmatched = true;

% Required
p.addOptional('sceneGeometryFileName', [], @(x)(isempty(x) || ischar(x)));

% parse
p.parse(sceneGeometryFileName)

% offer a file picker dialog if needed
if isempty(sceneGeometryFileName)
    [fileName, path] = uigetfile({'*_sceneGeometry.mat'});
    if isempty(fileName)
        returne
    end
    sceneGeometryFileName = [path, fileName];
end

% Load the sceneGeometry
dataLoad=load(sceneGeometryFileName);
sceneGeometry=dataLoad.sceneGeometry;
clear dataLoad

% Identify the frames of the ellipse array
ellipseArrayList = sceneGeometry.meta.estimateSceneParams.ellipseArrayList;
ellipseArrayList = ellipseArrayList(1);

% Set the initial scene params
sceneVector = sceneGeometry.meta.estimateSceneParams.search.x;

% Load the pupil perimeter data. It will be a structure variable
% "perimeter", with the fields .data and .meta
perimeterFileName = strrep(sceneGeometryFileName,'_sceneGeometry.mat','_correctedPerimeter.mat');
dataLoad=load(perimeterFileName);
perimeter=dataLoad.perimeter;
clear dataLoad

% Identify and open the corresponding video file
grayVideoName = strrep(sceneGeometryFileName,'_sceneGeometry.mat','_gray.avi');
videoInObj = videoIOWrapper(grayVideoName,'ioAction','read');

% Get the video properties
videoSizeX = videoInObj.Width;
videoSizeY = videoInObj.Height;
nFrames = floor(videoInObj.Duration*videoInObj.FrameRate);

% Pre-load the video frames we will use
% Define a variable to hold the selected frames
sourceFrames = zeros(videoSizeY,videoSizeX,3,length(ellipseArrayList),'uint8');
for ii = 1:length(ellipseArrayList)
    idx = ellipseArrayList(ii);
    % Obtain and render the frame
    videoInObj.CurrentTime = (idx - 1)/(videoInObj.FrameRate);
    sourceFrames(:,:,:,ii) = readFrame(videoInObj);
end

    % close the video object
    clear videoInObj

% Define a variable to hold the selected frames
framesToMontage = zeros(videoSizeY,videoSizeX,3,length(ellipseArrayList),'uint8');

% Turn off a warning that can occur during the montage step
warningState = warning;
warning('off','images:initSize:adjustingMag');

% Provide some instructions for the operator
fprintf('Adjust camera translation with i-j-k-l.\n');
fprintf('Adjust camera torsion with the left and right arrow keys.\n');
fprintf('Press return to be prompted to enter scene param values as text.\n');
fprintf('Press esc to exit.\n');

% Prepare a figure

% Enter a while loop
notDoneFlag = true;

while notDoneFlag
    for ii = 1:length(ellipseArrayList)
tic
        idx = ellipseArrayList(ii);
        
        % Obtain the eye pose from the boundary points from the perimeter
        Xp = perimeter.data{idx}.Xp;
        Yp = perimeter.data{idx}.Yp;
        eyePose = eyePoseEllipseFit(Xp, Yp, sceneGeometry, 'repeatSearchThresh', 10);
        
        % Add the rendered eye model
        if ~any(isnan(eyePose))
            sourceFrame = squeeze(sourceFrames(:,:,:,ii));
            [~,~,thisFrame] = renderEyePose(eyePose, sceneGeometry, ...
                'newFigure', true, 'visible', false, ...
                'backgroundImage',sourceFrame, ...
                'modelEyeLabelNames', {'retina' 'irisPerimeter' 'pupilEllipse' 'cornea'}, ...
                'modelEyePlotColors', {'.w' '.b' '-g' '.y'}, ...
                'modelEyeSymbolSizeScaler',1.5,...
                'modelEyeAlpha', 0.25);
        end
        % Add a text label for the frame number
        frameLabel = sprintf('frame: %d',idx);
        thisFrame.cdata = insertText(thisFrame.cdata,[20 20],frameLabel,'FontSize',30);
        % Store the frame
        framesToMontage(:,:,:,ii) = thisFrame.cdata;
    end

    
    % Prepare the main figure
    figHandle=figure('Visible','on');

    % Create the montage
    if length(ellipseArrayList)==1
        imshow(squeeze(framesToMontage(:,:,:,1)));
    else
        montage(framesToMontage);
    end
    drawnow
    toc
    pause
    
    close figHandle
end

end % GetWithDefault

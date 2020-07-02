function [ x, candidateSceneGeometry ] = estimateSceneParamsGUI(sceneGeometryFileName, varargin)
% Adjust scene parameter values
%
% Syntax:
%  [ initialParams ] = estimateSceneParamsGUI(grayVideoName)
%
% Description:
%
% Examples:
%{
    % ETTBSkip -- This is an idiosyncratic example.
    sceneGeometryFileName = '~/Dropbox (Aguirre-Brainard Lab)/TOME_processing/session1_restAndStructure/TOME_3045/042319/EyeTracking/GazeCal_sceneGeometry.mat';
    initialParams = estimateSceneParamsGUI(sceneGeometryFileName)
%}


%% Input parser
p = inputParser; p.KeepUnmatched = true;

% Required
p.addOptional('sceneGeometryFileName', [], @(x)(isempty(x) || ischar(x)));

% Optional
p.addParameter('frameSet',[],@isnumeric);

% parse
p.parse(sceneGeometryFileName,varargin{:})

% Load the sceneGeometry
dataLoad=load(sceneGeometryFileName);
sceneGeometry=dataLoad.sceneGeometry;
clear dataLoad

sceneGeometry.eye.landmarks.lateralCanthus.coords = [-10 -17 1];
sceneGeometry.eye.landmarks.medialCanthus.coords = [-3 12.5 -1];

% Identify the video file
grayVideoName = strrep(sceneGeometryFileName,'_sceneGeometry.mat','_gray.avi');

% Identify the frames of the ellipse array
if isempty(p.Results.frameSet)
    frameSet = sceneGeometry.meta.estimateSceneParams.obj.frameSet;
else
    frameSet = p.Results.frameSet;
end

% Get the video properties
videoInObj = videoIOWrapper(grayVideoName,'ioAction','read');
videoSizeX = videoInObj.Width;
videoSizeY = videoInObj.Height;

% Pre-load the video frames we will use
% Define a variable to hold the selected frames
sourceFrames = zeros(videoSizeY,videoSizeX,3,length(frameSet),'uint8');
for ii = 1:length(frameSet)
    frameIdx = frameSet(ii);
    % Obtain and render the frame
    videoInObj.CurrentTime = (frameIdx - 1)/(videoInObj.FrameRate);
    sourceFrames(:,:,:,ii) = readFrame(videoInObj);
end

% close the video object
clear videoInObj

% Load the pupil perimeter data. It will be a structure variable
% "perimeter", with the fields .data and .meta
perimeterFileName = strrep(sceneGeometryFileName,'_sceneGeometry.mat','_correctedPerimeter.mat');
dataLoad=load(perimeterFileName);
perimeter=dataLoad.perimeter;
clear dataLoad

% Load the glint data. It will be a structure variable
% "perimeter", with the fields .data and .meta
glintFileName = strrep(sceneGeometryFileName,'_sceneGeometry.mat','_glint.mat');
load(glintFileName,'glintData');

% Load the relativeCameraPosition file if it exists.
relativeCameraPositionFileName = strrep(sceneGeometryFileName,'_sceneGeometry.mat','_relativeCameraPosition.mat');
if exist(relativeCameraPositionFileName, 'file')==2
    dataLoad=load(relativeCameraPositionFileName);
    relativeCameraPosition=dataLoad.relativeCameraPosition;
    clear dataLoad
else
    relativeCameraPosition=[];
end


% Provide some instructions for the operator
fprintf('Adjust horizontal /vertical camera translation with the arrow keys.\n');
fprintf('Adjust depth camera translation with + and -.\n');
fprintf('Adjust camera torsion with j and k.\n');
fprintf('Adjust joint and differential rotation depth with g, b and f, h.\n');
fprintf('Move forward and backward in the ellipse frames with a and s\n')
fprintf('Press return to be prompted to enter scene param values as text.\n');
fprintf('Press esc to exit.\n');

% Set the current index and scene params
arrayIdx = 1;
x = [sceneGeometry.cameraPosition.torsion; ...
    sceneGeometry.cameraPosition.translation; ...
    sceneGeometry.eye.meta.rotationCenterScalers'];

% Prepare the main figure
figHandle=figure('Visible','on');
annotHandle=[];

% Enter a while loop
notDoneFlag = true;

while notDoneFlag
    
    % Update the scene parameters
    candidateSceneGeometry = sceneGeometry;
    candidateSceneGeometry.cameraPosition.torsion = x(1);
    candidateSceneGeometry.cameraPosition.translation = x(2:4);
    candidateSceneGeometry.eye.meta.rotationCenterScalers = x(5:6)';
    candidateSceneGeometry.eye.rotationCenters = human.rotationCenters( candidateSceneGeometry.eye );
    
    % Identify the frame to display
    frameIdx = frameSet(arrayIdx);
    
    % Adjust for relative camera position
    adjustedSceneGeometry = candidateSceneGeometry;
    if ~isempty(relativeCameraPosition)
        cameraPosition = candidateSceneGeometry.cameraPosition.translation;
        cameraPosition = cameraPosition - relativeCameraPosition.values(:,frameIdx);
        adjustedSceneGeometry.cameraPosition.translation = cameraPosition;
    end
    
    % Obtain the eye pose from the boundary points from the perimeter
    Xp = perimeter.data{frameIdx}.Xp;
    Yp = perimeter.data{frameIdx}.Yp;
    glintCoord = [glintData.X(frameIdx,:), glintData.Y(frameIdx,:)];
    
    eyePose = eyePoseEllipseFit(Xp, Yp, glintCoord, adjustedSceneGeometry,'cameraTransBounds',[0;0;0]);
    
    % Show this video frame
    thisFrame = sourceFrames(:,:,:,arrayIdx);
    frameLabel = sprintf('frame: %d',frameIdx);
    thisFrame = insertText(thisFrame,[20 20],frameLabel,'FontSize',30);
    imshow(squeeze(thisFrame));
    
    % Add the pupil perimeter points
    hold on
    plot(Xp ,Yp, '.w', 'MarkerSize', 1);
    plot(glintData.X(frameIdx,:), glintData.Y(frameIdx,:),'or');
    
    
    % Add the rendered eye model
    if ~any(isnan(eyePose))
        renderEyePose(eyePose, adjustedSceneGeometry, ...
            'newFigure', false, 'visible', true, ...
            'showAzimuthPlane', true, ...
            'modelEyeLabelNames', {'retina' 'irisPerimeter' 'pupilEllipse' 'cornea' 'aziRotationCenter' 'medialCanthus' 'lateralCanthus'}, ...
            'modelEyePlotColors', {'.w' '.b' '-g' '.y' '+c' 'Qg' 'Qg'}, ...
            'modelEyeSymbolSizeScaler',1.5,...
            'modelEyeAlpha', 0.25);
    end
    
    if ~isempty(annotHandle)
        delete(annotHandle)
    end
    
    % Wait for operator input
    fprintf('torsion: %0.2f, translation [%0.2f; %0.2f; %0.2f], rotation centers [%0.2f; %0.2f], eyePose [%0.2f; %0.2f; %0.2f; %0.2f]\n',x,eyePose);
    waitforbuttonpress
    keyChoiceValue = double(get(gcf,'CurrentCharacter'));
    switch keyChoiceValue
        case 28
            text_str = 'translate left';
            x(2)=x(2)-0.5;
        case 29
            text_str = 'translate right';
            x(2)=x(2)+0.5;
        case 30
            text_str = 'translate up';
            x(3)=x(3)+0.5;
        case 31
            text_str = 'translate down';
            x(3)=x(3)-0.5;
        case {45 95}
            text_str = 'translate closer in';
            x(4)=x(4)-1;
        case {61 43}
            text_str = 'translate farther away';
            x(4)=x(4)+1;
        case 103
            text_str = 'deeper rotation center';
            x(5)=x(5)+0.1;
        case 98
            text_str = 'shallower rotation center';
            x(5)=x(5)-0.1;
        case 102
            text_str = 'closer azi and ele rotation centers';
            x(6)=x(6)-0.1;
        case 104
            text_str = 'farther azi and ele rotation centers';
            x(6)=x(6)+0.1;
        case 97
            text_str = 'prior frame';
            arrayIdx = arrayIdx-1;
            if arrayIdx == 0
                arrayIdx = length(frameSet);
            end
        case 115
            text_str = 'next frame';
            arrayIdx = arrayIdx+1;
            if arrayIdx > length(frameSet)
                arrayIdx = 1;
            end
        case 106
            text_str = 'counter-clockwise torsion';
            x(1)=x(1)-1;
        case 107
            text_str = 'clockwise torsion';
            x(1)=x(1)+1;
        case 13
            x(1:4) = input('Enter values in square brackets, separated by semi-colons [tor;x;y;z]:');
            text_str = 'manual param entry';
        case 27
            notDoneFlag = false;
        otherwise
            text_str = 'unrecognized command';
    end
    
    if notDoneFlag
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
        hold off
    end
end

% Clean up
close(figHandle)
fprintf('\n');
fprintf('scene parameters = [%0.2f; %0.2f; %0.2f; %0.2f]\n',x);

end % Main function



%% LOCAL FUNCTION

function r = issubfield(s, f)
if isempty(f) || isempty(s)
    r = false;
else
    t = textscan(f,'%s','delimiter','.');
    t = t{1};
    r = true;
    for k = 1:numel(t)
        if isfield(s, t{k})
            s = s.(t{k});
        else
            r = false;
            return;
        end
    end
end
end % issubfield

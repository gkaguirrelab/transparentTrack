function makeEyeModelVideo(videoOutFileName,pupilFileName, sceneGeometryFileName, varargin)
% Create and store a video that displays the eye model fit to the data
%
% Syntax:
%  makeEyeModelVideo(videoOutFileName,pupilFileName, sceneGeometryFileName)
%
% Description:
%   This routine creates a video that illustrates for each frame the
%   appearance of the model eye in the image plane.
%
% Inputs:
%   videoOutFileName      - Full path to the output .avi file
%   pupilFileName         - Full path to a pupil data file. The file must
%                           have an eyePoses field.
%   sceneGeometryFileName - Full path to the sceneGeometry file
%
% Optional key/value pairs:
%  'verbose'              - Logical. Default false.
%  'videoOutFrameRate'    - Frame rate (in Hz) of saved video [default 60]
%  'saveCompressedVideo'  - Default value is true, resulting in a
%                           a video with a 10x reduction in file size
%  'modelEyeLabelNames'   - A cell array of the classes of eye model points
%                           to be displayed.
%  'modelEyePlotColors'   - The colors to be used for the plotting of each
%                           of the model eye label names.
%  'fitLabel'             - The field of the pupilData file to use
%
% Outputs:
%   None
%

%% Parse vargin for options passed here
p = inputParser; p.KeepUnmatched = true;

% Required
p.addRequired('videoOutFileName', @ischar);
p.addRequired('pupilFileName', @ischar);
p.addRequired('sceneGeometryFileName', @ischar);

% Optional display and I/O params
p.addParameter('verbose',false, @islogical);
p.addParameter('videoOutFrameRate', 60, @isnumeric);
p.addParameter('saveCompressedVideo', true, @islogical);
p.addParameter('modelEyeLabelNames', {'aziRotationCenter', 'eleRotationCenter', 'posteriorChamber' 'irisPerimeter' 'pupilPerimeter' 'anteriorChamber' 'cornealApex'}, @iscell);
p.addParameter('modelEyePlotColors', {'>r' '^m' '.w' '.b' '*g' '.y' '*y'}, @iscell);
p.addParameter('fitLabel', 'radiusSmoothed', @ischar);

% parse
p.parse(videoOutFileName, pupilFileName, sceneGeometryFileName, varargin{:})


%% Prepare variables

% Read in the pupilData file
dataLoad = load(p.Results.pupilFileName);
pupilData = dataLoad.pupilData;
clear dataLoad
eyePoses = pupilData.(p.Results.fitLabel).eyePoses.values;

% Load the sceneGeometry file
sceneGeometry = loadSceneGeometry(p.Results.sceneGeometryFileName, p.Results.verbose);

% Open a video object for writing
if p.Results.saveCompressedVideo
    videoOutObj = VideoWriter(videoOutFileName);
    videoOutObj.FrameRate = p.Results.videoOutFrameRate;
    open(videoOutObj);
else
    % Create a color map
    cmap = [linspace(0,1,256)' linspace(0,1,256)' linspace(0,1,256)'];
    cmap(1,:)=[1 0 0];
    cmap(2,:)=[0 1 0];
    cmap(3,:)=[0 0 1];
    cmap(4,:)=[1 1 0];
    cmap(5,:)=[0 1 1];
    cmap(6,:)=[1 0 1];
    
    videoOutObj = VideoWriter(videoOutFileName,'Indexed AVI');
    videoOutObj.FrameRate = p.Results.videoOutFrameRate;
    videoOutObj.Colormap = cmap;
    open(videoOutObj);
end

videoSizeX = sceneGeometry.cameraIntrinsic.sensorResolution(1);
videoSizeY = sceneGeometry.cameraIntrinsic.sensorResolution(2);

% A blank frame to initialize each frame
blankFrame = zeros(videoSizeY,videoSizeX)+0.5;

% Obtain the number of frames
nFrames = size(eyePoses,1);

% Open a figure
frameFig = figure( 'Visible', 'off');

% Alert the user
if p.Results.verbose
    tic
    fprintf(['Creating and saving model video. Started ' char(datetime('now')) '\n']);
    fprintf('| 0                      50                   100%% |\n');
    fprintf('.\n');
end


%% Loop through the frames
for ii = 1:nFrames
    
    % Update the progress display
    if p.Results.verbose && mod(ii,round(nFrames/50))==0
        fprintf('\b.\n');
    end
    
    % Plot the blank frame
    imshow(blankFrame, 'Border', 'tight');
    hold on
    axis off
    axis equal
    xlim([0 videoSizeX]);
    ylim([0 videoSizeY]);
    
    if ~any(isnan(eyePoses(ii,:)))
        
        % Obtain the pupilProjection of the model eye to the image plane
        [~, renderedFrame] = renderEyePose(eyePoses(ii,:), sceneGeometry, 'newFigure', false, 'modelEyeLabelNames', p.Results.modelEyeLabelNames, 'modelEyePlotColors', p.Results.modelEyePlotColors);
        
    end

    % Get ready for the next frame
    hold off
    
    % Write out this frame
    if p.Results.saveCompressedVideo
        thisFrame = squeeze(renderedFrame);
        writeVideo(videoOutObj,renderedFrame);
    else
        indexedFrame = rgb2ind(renderedFrame, cmap, 'nodither');
        writeVideo(videoOutObj,indexedFrame);
    end
    
end % Loop over frames


%% Save and cleanup

% Close the figure
close(frameFig);

% close the video objects
clear videoOutObj videoInObj

% report completion of fit video generation
if p.Results.verbose
    toc
    fprintf('\n');
end


end % function

function makeEyeModelVideo(videoOutFileName,pupilFileName, sceneGeometryFileName, varargin)
% Create and store a video that displays the eye model fit to the data
%
% Description:
%   This routine creates a video that illustrates for each frame the
%   appearance of the model eye in the image plane.
%
% Inputs:
%   videoOutFileName      - Full path to the output .avi file
%   pupilFileName         - Full path to a pupil data file. The file must
%                           have an eyeParams field.
%   sceneGeometryFileName - Full path to the sceneGeometry file
%
% Optional key/value pairs (display and I/O):
%  'verbosity'            - Level of verbosity. [none, full]
%  'videoOutFrameRate'    - Frame rate (in Hz) of saved video [default 60]
%  'saveCompressedVideo'  - Default value is true, resulting in a
%                           a video with a 10x reduction in file size
%  'videoSizeX'           - Size of the video in the X dimension
%  'videoSizeY'           - Size of the video in the Y dimension
%  'labelNames'           - A cell array of the classes of eye model points
%                           to be displayed.
%  'plotColors'           - The colors to be used for the plotting of each
%                           of the label names.
%  'ellipseFitLabel'      - The field of the pupilData file that contains
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
p.addParameter('verbosity','none', @ischar);
p.addParameter('videoOutFrameRate', 60, @isnumeric);
p.addParameter('saveCompressedVideo', true, @islogical);
p.addParameter('videoSizeX', 640, @isnumeric);
p.addParameter('videoSizeY', 480, @isnumeric);
p.addParameter('labelNames', {'rotationCenter', 'posteriorChamber' 'irisPerimeter' 'pupilPerimeter' 'anteriorChamber'}, @iscell);
p.addParameter('plotColors', {'+r' '.w' '.b' '.g' '.y'}, @iscell);
p.addParameter('ellipseFitLabel', 'radiusSmoothed', @ischar);

% parse
p.parse(videoOutFileName, pupilFileName, sceneGeometryFileName, varargin{:})


%% Alert the user and prepare variables
if strcmp(p.Results.verbosity,'full')
    tic
    fprintf(['Creating and saving model video. Started ' char(datetime('now')) '\n']);
    fprintf('| 0                      50                   100%% |\n');
    fprintf('.\n');
end

% Read in the pupilData file
dataLoad = load(p.Results.pupilFileName);
pupilData = dataLoad.pupilData;
clear dataLoad
eyeParams = pupilData.(p.Results.ellipseFitLabel).eyeParams.values;

% Read in the sceneGeometry file
dataLoad = load(p.Results.sceneGeometryFileName);
sceneGeometry = dataLoad.sceneGeometry;
clear dataLoad

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

% A blank frame to initialize each frame
blankFrame = zeros(480,640)+0.5;

% Obtain the number of frames
nFrames = size(eyeParams,1);

% Open a figure
frameFig = figure( 'Visible', 'off');

%% Loop through the frames
for ii = 1:nFrames
    
    % Update the progress display
    if strcmp(p.Results.verbosity,'full') && mod(ii,round(nFrames/50))==0
        fprintf('\b.\n');
    end
    
    % Plot the blank frame
    imshow(blankFrame, 'Border', 'tight');
    hold on
    axis off
    axis equal
    xlim([0 p.Results.videoSizeX]);
    ylim([0 p.Results.videoSizeY]);
    
    if ~any(isnan(eyeParams(ii,:)))
        
        % Obtain the pupilProjection of the model eye to the image plane
        [~, ~, imagePoints, pointLabels] = pupilProjection_fwd(eyeParams(ii,:), sceneGeometry, true);
        
        % Loop through the point labels present in the eye model
        for pp = 1:length(p.Results.labelNames)
            idx = strcmp(pointLabels,p.Results.labelNames{pp});
            plot(imagePoints(idx,1), imagePoints(idx,2), p.Results.plotColors{pp})
        end
        
    end
    
    % Clean up the plot
    hold off
    
    % Get the frame and close the figure
    thisFrame=getframe(frameFig);
    
    % Write out this frame
    if p.Results.saveCompressedVideo
        thisFrame = squeeze(thisFrame);
        writeVideo(videoOutObj,thisFrame);
    else
        indexedFrame = rgb2ind(thisFrame, cmap, 'nodither');
        writeVideo(videoOutObj,indexedFrame);
    end
    
end % Loop over frames


%% Save and cleanup

% Close the figure
close(frameFig);

% close the video objects
clear videoOutObj videoInObj

% report completion of fit video generation
if strcmp(p.Results.verbosity,'full')
    toc
    fprintf('\n');
end


end % function

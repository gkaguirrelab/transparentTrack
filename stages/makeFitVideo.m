function makeFitVideo(videoInFileName, videoOutFileName, varargin)
% Create and store a video that displays the results of eye tracking
%
% Syntax:
%  makeFitVideo(videoInFileName, videoOutFileName)
%
% Description:
%   This routine creates an integrated fit video that illustrates the
%   position of the pupil and glint and indicates any control instructions.
%
% Inputs:
%	videoInFileName       - Full path to an .avi file; typically the "gray"
%                           video
%   videoOutFileName      - Full path to the output .avi file
%
% Optional key/value pairs (display and I/O):
%  'verbose'              - Logical. Default false.
%  'videoOutFrameRate'    - Frame rate (in Hz) of saved video [default 60]
%  'saveCompressedVideo'  - Default value is true, resulting in a
%                           a video with a 10x reduction in file size
%  'nFrames'              - Analyze fewer than the total number of frames.
%  'glint/perimeter/pupil/sceneGeometry/FileName' - Full path to a file
%                           to be included in the video. 
%  'glint/perimeter/pupil/sceneGeometry/Color' - Text string that assigns
%                           a color to the display of this item.
%  'fitLabel'             - The field of the pupilData file that contains
%                           ellipse fit params to be added to the video.
%  'controlFileName'      - Full path to the control file to be included.
%
% Outputs:
%   None
%

%% Parse vargin for options passed here
p = inputParser; p.KeepUnmatched = true;

% Required
p.addRequired('videoInFileName', @ischar);
p.addRequired('videoOutFileName', @ischar);

% Optional display and I/O params
p.addParameter('verbose',false,@islogical);
p.addParameter('videoOutFrameRate', 60, @isnumeric);
p.addParameter('saveCompressedVideo', true, @islogical);

% Optional flow control params
p.addParameter('nFrames',Inf,@isnumeric);

% Optional video items
p.addParameter('glintFileName',[],@(x)(isempty(x) | ischar(x)));
p.addParameter('perimeterFileName',[],@(x)(isempty(x) | ischar(x)));
p.addParameter('pupilFileName',[],@(x)(isempty(x) | ischar(x)));
p.addParameter('sceneGeometryFileName',[],@(x)(isempty(x) | ischar(x)));
p.addParameter('glintColor','r',@ischar);
p.addParameter('perimeterColor','w',@ischar);
p.addParameter('pupilColor','green',@ischar);
p.addParameter('sceneGeometryColor','magenta',@ischar);
p.addParameter('modelEyeLabelNames', {'posteriorChamber' 'irisPerimeter' 'anteriorChamber'}, @iscell);
p.addParameter('modelEyePlotColors', {'.w' '.b' '.y'}, @iscell);
p.addParameter('modelEyeAlpha', 0,@isnumeric);
p.addParameter('fitLabel', 'radiusSmoothed',@(x)(isempty(x) | ischar(x)));
p.addParameter('controlFileName',[],@(x)(isempty(x) | ischar(x)));

% parse
p.parse(videoInFileName, videoOutFileName, varargin{:})


%% Prepare variables

% Read in the glint file if passed
if ~isempty(p.Results.glintFileName)
    dataLoad = load(p.Results.glintFileName);
    glintData = dataLoad.glintData;
    clear dataLoad
else
    glintData=[];
end

% Read in the perimeter file if passed
if ~isempty(p.Results.perimeterFileName)
    dataLoad = load(p.Results.perimeterFileName);
    perimeter = dataLoad.perimeter;
    clear dataLoad
else
    perimeter=[];
end

% Read in the pupilData file if passed
if ~isempty(p.Results.pupilFileName)
    dataLoad = load(p.Results.pupilFileName);
    pupilData = dataLoad.pupilData;
    clear dataLoad
    ellipseFitParams = pupilData.(p.Results.fitLabel).ellipses.values;
    if isfield(pupilData.(p.Results.fitLabel),'eyePoses')
        eyePoses = pupilData.(p.Results.fitLabel).eyePoses.values;
    else
        eyePoses = [];
    end
else
    ellipseFitParams=[];
    eyePoses=[];
end

% Read in and parse the control file if passed
if ~isempty(p.Results.controlFileName)
    instructions = loadControlFile(p.Results.controlFileName);
else
    instructions(1).frame=[];
    instructions(1).type=[];
    instructions(1).params=[];
end

% Load the sceneGeometry file if passed
if ~isempty(p.Results.sceneGeometryFileName)
    sceneGeometry = loadSceneGeometry(p.Results.sceneGeometryFileName, p.Results.verbose);
else
    sceneGeometry = [];
end

% Open a video object for reading
videoInObj = VideoReader(videoInFileName);

% get number of frames
if p.Results.nFrames == Inf
    nFrames = floor(videoInObj.Duration*videoInObj.FrameRate);
else
    nFrames = p.Results.nFrames;
end

% get video dimensions
videoSizeX = videoInObj.Width;
videoSizeY = videoInObj.Height;

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

% Create blank figure. Save the the axis handles. By allowing the figure to
% be visible, processing is about 10x faster. This is a weird property of
% how MATLAB handles frame buffering. 
hFig = figure( 'Visible', 'on');
hAxes = gca();

% Alert the user
if p.Results.verbose
    tic
    fprintf(['Creating and saving fit video. Started ' char(datetime('now')) '\n']);
    fprintf('| 0                      50                   100%% |\n');
    fprintf('.\n');
end


%% Loop through the frames
for ii = 1:nFrames
    
    % Update the progress display
    if p.Results.verbose && mod(ii,round(nFrames/50))==0
        fprintf('\b.\n');
    end
    
    % read the source video frame into memory
    sourceFrame = readFrame(videoInObj);
    sourceFrame = rgb2gray (sourceFrame);
            
    % show the initial frame
    if ii==1
        hImage = imshow(sourceFrame,'Border', 'tight','Parent',hAxes);
    else
        hold off
        set(hImage,'CData',sourceFrame);
    end
    hold on
    
    hPlot = gobjects(0);
    
    % add glint
    if ~isempty(p.Results.glintFileName)
        for gg = 1:size(glintData.X,2)
        hPlot(end+1) = plot(glintData.X(ii,gg),glintData.Y(ii,gg),['*' p.Results.glintColor]);
        end
    end
    
    % add pupil perimeter
    if ~isempty(p.Results.perimeterFileName)
        % get the data frame
        if ~isempty(perimeter.data{ii}.Xp)
            hPlot(end+1) = plot(perimeter.data{ii}.Xp ,perimeter.data{ii}.Yp, ['.' p.Results.perimeterColor], 'MarkerSize', 1);
        end
    end
        
    % add pupil ellipse fit
    if ~isempty(p.Results.pupilFileName)
        if ~isempty(ellipseFitParams)
            if sum(isnan(ellipseFitParams(ii,:)))==0
                % build ellipse impicit equation
                pFitImplicit = ellipse_ex2im(ellipse_transparent2ex(ellipseFitParams(ii,:)));
                fh=@(x,y) pFitImplicit(1).*x.^2 +pFitImplicit(2).*x.*y +pFitImplicit(3).*y.^2 +pFitImplicit(4).*x +pFitImplicit(5).*y +pFitImplicit(6);
                % superimpose the ellipse using fimplicit or ezplot (ezplot
                % is the fallback option for older Matlab versions)
                if exist('fimplicit','file')==2
                    hPlot(end+1) = fimplicit(fh,[1, videoSizeX, 1, videoSizeY],'Color', p.Results.pupilColor,'LineWidth',1);
                    set(gca,'position',[0 0 1 1],'units','normalized')
                    axis off;
                else
                    hPlot(end+1) = ezplot(fh,[1, videoSizeX, 1, videoSizeY]);
                    set(hPlot(end), 'Color', p.Results.pupilColor)
                    set(hPlot(end), 'LineWidth',1);
                end
            end
        end
    end

    % superimpose the model eye
    if ~isempty(eyePoses) && sum(p.Results.modelEyeAlpha)~=0
        % If we have a defined eyePose for this frame, display the modelEye
        if ~any(isnan(eyePoses(ii,:)))
            [~, hRender] = renderEyePose(eyePoses(ii,:), sceneGeometry, 'newFigure', false, ...
                'modelEyeLabelNames', p.Results.modelEyeLabelNames, ...
                'modelEyePlotColors', p.Results.modelEyePlotColors, ...
                'modelEyeAlpha', p.Results.modelEyeAlpha);   
        end
    end

    % add an instruction label
    if ~isempty(p.Results.controlFileName)
        instructionIdx = find ([instructions.frame] == ii);
        if ~isempty(instructionIdx)
            text_str = instructions(instructionIdx(end)).type;
            hPlot(end+1) = annotation('textbox',...
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
        end
    end
        
    % Get the frame
    drawnow;
    thisFrame=getframe(hFig);
    
    % Write out this frame
    if p.Results.saveCompressedVideo
        thisFrame = squeeze(thisFrame);
        writeVideo(videoOutObj,thisFrame);
    else
        indexedFrame = rgb2ind(thisFrame, cmap, 'nodither');
        writeVideo(videoOutObj,indexedFrame);
    end
    
    % Clear the plot objects
    if exist('hRender', 'var')
        delete(hRender);
    end
    delete(hPlot);

        
end % Loop over frames


%% Save and cleanup

% Close the figure
close(hFig);

% close the video objects
clear videoOutObj videoInObj

% report completion of fit video generation
if p.Results.verbose
    toc
    fprintf('\n');
end


end % function

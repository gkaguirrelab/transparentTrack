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
%
% Optional key/value pairs (flow control):
%  'nFrames'              - Analyze fewer than the total number of frames.
%
% Optional key/value pairs (video items):
%  'glint/perimeter/pupil/sceneGeometry/controlFileName' - Full path to a 
%                           file to be included in the video. 
%  'adjustedCameraPositionTranslation' - 3x1 vector that provides position
%                           of the camera relative to the origin of the
%                           world coordinate system (which is the anterior
%                           surface of the cornea in primary gaze). This
%                           value is used to update the sceneGeometry file
%                           to account for head movement that has taken
%                           place between the sceneGeometry acquisition and
%                           the acquisition undergoing analysis. This
%                           updated camera position should reflect the
%                           camera position at the start of the current
%                           acquisition.
%  'relativeCameraPositionFileName' - Char. This is the full path to a
%                           relativeCameraPosition.mat file that provides
%                           the relative position of the camera at each
%                           video frame relative to the initial position of
%                           the camera.
%  'glint/perimeter/pupil/sceneGeometry/Color' - Text string that assigns
%                           a color to the display of this item.
%  'pupilRMSERangeAlphaScaler' - 1x2 vector. The alpha transparency of the
%                           fitted pupil ellipse varies between 0.25
%                           (somewhat transparent) to 1 (solid). The
%                           alpha level is determined by the RMSE of the
%                           ellipse fit to the pupil perimeter, with the
%                           mapping determined by the minimum and maximum
%                           RMSE values provided in this vector.
%  'modelEyeLabelNames/modelEyePlotColors' - Cell arrays that are passed to
%                           renderEyePose to determine the model eye
%                           elements that are superimposed on the video.
%  'modelEyeMaxAlpha'     - Scalar. The alpha value of the superimposed
%                           model eye. Defaults to zero, which results in
%                           the model eye not being displayed.
%  'modelEyeRMSERangeAlphaScaler' - 1x2 vector. Maps the RMSE fit of the
%                           pupil ellipse to the alpha value of the model
%                           eye.
%  'modelEyeSymbolSizeScaler' - Scalar. Determines the size of the plot
%                           symbols used to render the model eye.
%  'suppressBlinks'       - Logical. If set to true, the model eye will not
%                           be rendered for frames marked as "blinks" in 
%                           the instructions.
%  'uniformityThreshold'  - Scalar. A uniformity score is calculated for
%                           each frame that describes the extent to which
%                           the pupil perimeter points are either minimally
%                           (0) or fully (1) distributed around the
%                           ellipse. The eyeModel is only displayed for
%                           those frames in which the linear uniformity
%                           value is above this threshold.
%  'fitLabel'             - The field of the pupilData file that contains
%                           ellipse fit params to be added to the video.
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
p.addParameter('verbose', false, @islogical);
p.addParameter('videoOutFrameRate', 60, @isnumeric);
p.addParameter('saveCompressedVideo', true, @islogical);

% Optional flow control params
p.addParameter('nFrames', Inf, @isnumeric);

% Optional video items
p.addParameter('glintFileName', [], @(x)(isempty(x) | ischar(x)));
p.addParameter('perimeterFileName', [], @(x)(isempty(x) | ischar(x)));
p.addParameter('pupilFileName', [], @(x)(isempty(x) | ischar(x)));
p.addParameter('sceneGeometryFileName', [], @(x)(isempty(x) | ischar(x)));
p.addParameter('controlFileName', [], @(x)(isempty(x) | ischar(x)));
p.addParameter('adjustedCameraPositionTranslation',[],@isnumeric);
p.addParameter('relativeCameraPositionFileName',[],@ischar);
p.addParameter('glintColor', 'r', @ischar);
p.addParameter('perimeterColor', 'w', @ischar);
p.addParameter('pupilColor', 'green', @ischar);
p.addParameter('pupilRMSERangeAlphaScaler',[1,4],@isnumeric);
p.addParameter('modelEyeLabelNames', {'retina' 'irisPerimeter' 'cornea' 'glint'}, @iscell);
p.addParameter('modelEyePlotColors', {'.w' '.b' '.y' 'or'}, @iscell);
p.addParameter('modelEyeMaxAlpha', 0, @isnumeric);
p.addParameter('modelEyeRMSERangeAlphaScaler',[1,4],@isnumeric);
p.addParameter('modelEyeSymbolSizeScaler',1,@isnumeric);
p.addParameter('suppressBlinks',true,@islogical);
p.addParameter('uniformityThreshold',0.33,@isscalar);
p.addParameter('fitLabel', 'radiusSmoothed', @(x)(isempty(x) | ischar(x)));

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
ellipseFitParams=[];
ellipseFitRMSE=[];
eyePoses=[];
linearUniformity = [];
fitAtBound = [];
if ~isempty(p.Results.pupilFileName)
    dataLoad = load(p.Results.pupilFileName);
    pupilData = dataLoad.pupilData;
    clear dataLoad
    ellipseFitParams = pupilData.(p.Results.fitLabel).ellipses.values;
    ellipseFitRMSE = pupilData.(p.Results.fitLabel).ellipses.RMSE;
    if isfield(pupilData.(p.Results.fitLabel),'eyePoses')
        eyePoses = pupilData.(p.Results.fitLabel).eyePoses.values;
    end
    % Get the uniformity field if is available
    if isfield(pupilData.(p.Results.fitLabel).ellipses,'uniformity')
        linearUniformity = pupilData.(p.Results.fitLabel).ellipses.uniformity;
    end
    % Get the fitAtBound field. Use the eyePose vector if available
    fitAtBound = [];
    if isfield(pupilData.(p.Results.fitLabel).ellipses,'fitAtBound')
        fitAtBound = pupilData.(p.Results.fitLabel).ellipses.fitAtBound;
    end
    if isfield(pupilData.(p.Results.fitLabel),'eyePoses')
        if isfield(pupilData.(p.Results.fitLabel).eyePoses,'fitAtBound')
            fitAtBound = pupilData.(p.Results.fitLabel).eyePoses.fitAtBound;
        end
    end
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
    % Load the sceneGeometry file
    dataLoad=load(p.Results.sceneGeometryFileName);
    sceneGeometry=dataLoad.sceneGeometry;
    clear dataLoad
    % An earlier version of the code defined a non-zero iris thickness. We
    % force this to zero here to speed computation
    sceneGeometry.eye.iris.thickness=0;
    % If an adjustedCameraPositionTranslation value has been passed, update this field
    % of the sceneGeometry
    if ~isempty(p.Results.adjustedCameraPositionTranslation)
        sceneGeometry.cameraPosition.translation = p.Results.adjustedCameraPositionTranslation;
    end
else
    sceneGeometry = [];
end

% Load the relativeCameraPosition file if passed and it exists
if ~isempty(p.Results.relativeCameraPositionFileName)
    if exist(p.Results.relativeCameraPositionFileName, 'file')==2
        dataLoad=load(p.Results.relativeCameraPositionFileName);
        relativeCameraPosition=dataLoad.relativeCameraPosition;
        clear dataLoad
    else
        relativeCameraPosition=[];
    end
else
    relativeCameraPosition=[];
end

% Prepare the video object
videoInObj = videoIOWrapper(videoInFileName,'ioAction','read');

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
    videoOutObj = videoIOWrapper(videoOutFileName,'ioAction','write');
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
    
    videoOutObj = videoIOWrapper(videoOutFileName,'ioAction','write','indexedAVI',true);
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
        hPlot(end+1) = plot(glintData.X(ii,gg),glintData.Y(ii,gg),['.' p.Results.glintColor]);
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
                % If the pupilEllipse (or corresponding eyePose ellipse)
                % hit a bound, then force the pupil plot color to be red
                pupilColor = p.Results.pupilColor;
                if ~isempty(fitAtBound)
                    if fitAtBound(ii)
                        pupilColor = 'red';
                    end
                end
                % build ellipse impicit equation
                pFitImplicit = ellipse_ex2im(ellipse_transparent2ex(ellipseFitParams(ii,:)));
                fh=@(x,y) pFitImplicit(1).*x.^2 +pFitImplicit(2).*x.*y +pFitImplicit(3).*y.^2 +pFitImplicit(4).*x +pFitImplicit(5).*y +pFitImplicit(6);
                % superimpose the ellipse using fimplicit or ezplot (ezplot
                % is the fallback option for older Matlab versions)
                if exist('fimplicit','file')==2
                    hPlot(end+1) = fimplicit(fh,[1, videoSizeX, 1, videoSizeY],'Color', pupilColor,'LineWidth',1);
                    set(gca,'position',[0 0 1 1],'units','normalized')
                    axis off;
                else
                    hPlot(end+1) = ezplot(fh,[1, videoSizeX, 1, videoSizeY]);
                    set(hPlot(end), 'Color', pupilColor)
                    set(hPlot(end), 'LineWidth',1);
                end
                % To support alpha transparency, we need to re-plot the
                % ellipse as a conventional MATLAB line object, and not as
                % an implicit function or a contour.
                
                xData = hPlot(end).XData;
                yData = hPlot(end).YData;
                delete(hPlot(end));
                if ~isempty(xData)
                    hPlot(end)=plot(xData,yData,'Color', pupilColor,'LineWidth',1);
                    % Scale the ellipse line alpha by the RMSE ellipse fit
                    % value for this frame
                    RMSEVal = max([ellipseFitRMSE(ii) p.Results.modelEyeRMSERangeAlphaScaler(1)]);
                    RMSEVal = min([RMSEVal p.Results.modelEyeRMSERangeAlphaScaler(2)]);
                    alphaVal = 1 - ( (RMSEVal-p.Results.pupilRMSERangeAlphaScaler(1))/(p.Results.pupilRMSERangeAlphaScaler(2)-p.Results.pupilRMSERangeAlphaScaler(1)));
                    alphaVal = max([alphaVal 0.25]);
                    hPlot(end).Color(4) = alphaVal;
                end
            end
        end
    end

    %% Superimpose the model eye
    % If an instruction field is available, and we should suppress blinks,
    % check if this frame is a blink
    notBlinkFrame = true;
    if ~isempty(p.Results.pupilFileName)
        if isfield(pupilData,'instructions') && p.Results.suppressBlinks
            notBlinkFrame = ~pupilData.instructions.blink(ii);
        end
    end
    % If a linear uniformity score is available for the frame, check to see
    % that the threshold is met to display the eye model
    sufficientlyUniform = true;
    if ~isempty(linearUniformity)
        sufficientlyUniform = linearUniformity(ii) > p.Results.uniformityThreshold;
    end
    if ~isempty(eyePoses) && sum(p.Results.modelEyeMaxAlpha)~=0 && notBlinkFrame && sufficientlyUniform
        % If a relativeCameraPosition is defined, update the
        % sceneGeometry
        adjustedSceneGeometry = sceneGeometry;
        if ~isempty(relativeCameraPosition)
            cameraPosition = sceneGeometry.cameraPosition.translation;
            cameraPosition = cameraPosition + relativeCameraPosition.values(:,ii);
            adjustedSceneGeometry.cameraPosition.translation = cameraPosition;
        end
        % Scale the model eye alpha by the RMSE ellipse fit value for this
        % frame
        RMSEVal = max([ellipseFitRMSE(ii) p.Results.modelEyeRMSERangeAlphaScaler(1)]);
        RMSEVal = min([RMSEVal p.Results.modelEyeRMSERangeAlphaScaler(2)]);
        alphaVal = p.Results.modelEyeMaxAlpha - p.Results.modelEyeMaxAlpha*( (RMSEVal-p.Results.modelEyeRMSERangeAlphaScaler(1))/(p.Results.modelEyeRMSERangeAlphaScaler(2)-p.Results.modelEyeRMSERangeAlphaScaler(1)));
        % If we have a defined eyePose for this frame, display the modelEye
        if ~any(isnan(eyePoses(ii,:)))
            [~, hRender] = renderEyePose(eyePoses(ii,:), adjustedSceneGeometry, 'newFigure', false, ...
                'modelEyeLabelNames', p.Results.modelEyeLabelNames, ...
                'modelEyePlotColors', p.Results.modelEyePlotColors, ...
                'modelEyeAlpha', alphaVal, ...
                'modelEyeSymbolSizeScaler', p.Results.modelEyeSymbolSizeScaler);
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
    thisFrame=getframe(hAxes);
    
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

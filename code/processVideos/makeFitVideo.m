function makeFitVideo(videoInFileName, videoOutFileName, varargin)
% makeFitVideo(perimeterVideoFileName, varargin)
%
% This routine creates an integrated fit video that illustrates the
% position of the pupil and glint and indicates any control instructions.
%
% INPUTS:
%   videoInFileName: full path to an .avi file. Typically the gray video
%   videoOutFileName: full path to an .avi file. This will be the output
%
% Optional key/value pairs (display and I/O)
%  'verbosity' - level of verbosity. [none, full]
%  'videoOutFrameRate' - frame rate (in Hz) of saved video
%
% Optional key/value pairs (flow control)
%  'nFrames' - analyze fewer than the total number of frames.
%  'useParallel' - If set to true, use the Matlab parallel pool for the
%    initial ellipse fitting.
%  'nWorkers' - Specify the number of workers in the parallel pool. If
%    undefined the default number will be used.
%  'tbtbProjectName' - The workers in the parallel pool are configured by
%    issuing a tbUseProject command for the project specified here.
%
% Optional items to include in the video
%	glintFileName
%   glintColor - indicate color with a single character or string (e.g.,
%       'r' for red or 'green')
%   perimeterFileName
%   perimeterColor
%   pupilFileName
%   pupilColor
%   whichFieldToPlot - The name of the field within the ellipseFitData
%   	struct that is to be plotted
%   irisFileName
%   irisColor
%   controlFileName
%

%% Parse vargin for options passed here
p = inputParser; p.KeepUnmatched = true;

% Required
p.addRequired('videoInFileName', @ischar);
p.addRequired('videoOutFileName', @ischar);

% Optional display and I/O params
p.addParameter('verbosity','none', @ischar);
p.addParameter('videoOutFrameRate', 60, @isnumeric);

% Optional flow control params
p.addParameter('nFrames',Inf,@isnumeric);
p.addParameter('useParallel',false,@islogical);
p.addParameter('nWorkers',[],@(x)(isempty(x) | isnumeric(x)));
p.addParameter('tbtbRepoName','transparentTrack',@ischar);

% Optional items to include in the video
p.addParameter('glintFileName',[],@(x)(isempty(x) | ischar(x)));
p.addParameter('glintColor','r',@ischar);
p.addParameter('perimeterFileName',[],@(x)(isempty(x) | ischar(x)));
p.addParameter('perimeterColor','w',@ischar);
p.addParameter('pupilFileName',[],@(x)(isempty(x) | ischar(x)));
p.addParameter('pupilColor','green',@ischar);
p.addParameter('whichFieldToPlot', 'ellipseParamsAreaSmoothed_mean',@(x)(isempty(x) | ischar(x)));
p.addParameter('irisFileName',[],@(x)(isempty(x) | ischar(x)));
p.addParameter('irisColor','red',@ischar);
p.addParameter('controlFileName',[],@(x)(isempty(x) | ischar(x)));
p.addParameter('sceneGeometryFileName',[],@(x)(isempty(x) | ischar(x)));
p.addParameter('sceneGeometryColor','magenta',@ischar);


% parse
p.parse(videoInFileName, videoOutFileName, varargin{:})


%% Set up the parallel pool
if p.Results.useParallel
    if strcmp(p.Results.verbosity,'full')
        tic
        fprintf(['Opening parallel pool. Started ' char(datetime('now')) '\n']);
    end
    if isempty(p.Results.nWorkers)
        parpool;
    else
        parpool(p.Results.nWorkers);
    end
    poolObj = gcp;
    if isempty(poolObj)
        nWorkers=0;
    else
        nWorkers = poolObj.NumWorkers;
        % Use TbTb to configure the workers.
        if ~isempty(p.Results.tbtbRepoName)
            spmd
                tbUse(p.Results.tbtbRepoName,'reset','full','verbose',false,'online',false);
            end
            if strcmp(p.Results.verbosity,'full')
                fprintf('CAUTION: Any TbTb messages from the workers will not be shown.\n');
            end
        end
    end
    if strcmp(p.Results.verbosity,'full')
        toc
        fprintf('\n');
    end
else
    nWorkers=0;
end


%% Alert the user and prepare variables
if strcmp(p.Results.verbosity,'full')
    tic
    fprintf(['Creating and saving fit video. Started ' char(datetime('now')) '\n']);
    fprintf('| 0                      50                   100%% |\n');
    fprintf('.\n');
end

% Read in the variables to display. If not loaded, explicitly set the
% variable to empty so that the parfor loop does not panic about uncalled
% lines of code that make reference to non-existent variables.

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
    pupilFitParams = pupilData.(p.Results.whichFieldToPlot);
else
    pupilFitParams=[];
end

% Read in the irisData file if passed
if ~isempty(p.Results.irisFileName)
    dataLoad = load(p.Results.irisFileName);
    irisData = dataLoad.irisData;
    clear dataLoad
    irisFitParams = irisData.pEllipseFitTransparent;
else
    irisFitParams=[];
end

% Read in and parse the control file if passed
if ~isempty(p.Results.controlFileName)
    instructions = loadControlFile(p.Results.controlFileName);
else
    instructions(1).frame=[];
    instructions(1).type=[];
    instructions(1).params=[];
end

% Read in the sceneGeometry file if passed
if ~isempty(p.Results.sceneGeometryFileName)
    dataLoad = load(p.Results.sceneGeometryFileName);
    sceneGeometry = dataLoad.sceneGeometry;
    clear dataLoad
else
    sceneGeometry=[];
end

% read video file into memory
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
% initialize variable to hold the perimeter data
sourceVideo = zeros(videoSizeY,videoSizeX,nFrames,'uint8');
% read the video into memory
for ii = 1:nFrames
    thisFrame = readFrame(videoInObj);
    sourceVideo(:,:,ii) = rgb2gray (thisFrame);
end
% close the video object
clear videoInObj

% create a temp directory for the output frames
scratchDirectoryName = [videoOutFileName '_scratch'];
if ~exist(scratchDirectoryName,'dir')
    mkdir(scratchDirectoryName)
end

% get glintData ready for the parfor. This includes transposing the
% variables
if ~isempty(p.Results.glintFileName)
    glintData_X = glintData.X;
    glintData_Y = glintData.Y;
else
    glintData_X = nan(1,nFrames);
    glintData_Y = nan(1,nFrames);
end

% Recast perimeter and the video source into a sliced cell array to reduce par for
% broadcast overhead
if ~isempty(perimeter)
    frameCellArray = perimeter.data(1:nFrames);
    clear perimeter
end

sourceVideoArray = arrayfun(@(ii) {squeeze(sourceVideo(:,:,ii))},1:1:nFrames);

%% Loop through the frames
parfor (ii = 1:nFrames, nWorkers)
    
    % Update the progress display
    if strcmp(p.Results.verbosity,'full') && mod(ii,round(nFrames/50))==0
        fprintf('\b.\n');
    end
    
    % Create a figure
    frameFig = figure( 'Visible', 'off');
    
    % show the initial frame
    imshow(sourceVideoArray{ii}, 'Border', 'tight');
    hold on
    
    % add glint
    if ~isempty(p.Results.glintFileName)
        plot(glintData_X(ii),glintData_Y(ii),['*' p.Results.glintColor]);
    end
    
    % add pupil perimeter
    if ~isempty(p.Results.perimeterFileName)
        % get the data frame
        if ~isempty(frameCellArray{ii}.Xp)
            plot(frameCellArray{ii}.Xp,frameCellArray{ii}.Yp,['.' p.Results.perimeterColor], 'MarkerSize', 1);
        end
    end
    
    % add pupil ellipse fit
    if ~isempty(p.Results.pupilFileName)
        if ~isempty(pupilFitParams)
            if sum(isnan(pupilFitParams(ii,:)))==0
                % build ellipse impicit equation
                pFitImplicit = ellipse_ex2im(ellipse_transparent2ex(pupilFitParams(ii,:)));
                fh=@(x,y) pFitImplicit(1).*x.^2 +pFitImplicit(2).*x.*y +pFitImplicit(3).*y.^2 +pFitImplicit(4).*x +pFitImplicit(5).*y +pFitImplicit(6);
                % superimpose the ellipse using fimplicit or ezplot (ezplot
                % is the fallback option for older Matlab versions)
                if exist('fimplicit','file')==2
                    fimplicit(fh,[1, videoSizeX, 1, videoSizeY],'Color', p.Results.pupilColor,'LineWidth',1);
                    set(gca,'position',[0 0 1 1],'units','normalized')
                    axis off;
                else
                    plotHandle=ezplot(fh,[1, videoSizeX, 1, videoSizeY]);
                    set(plotHandle, 'Color', p.Results.pupilColor)
                    set(plotHandle,'LineWidth',1);
                end
            end
        end
    end
    
    % add iris ellipse fit
    if ~isempty(p.Results.irisFileName)
        if ~isempty(irisFitParams)
            if sum(isnan(irisFitParams(ii,:)))==0
                % build ellipse impicit equation
                pFitImplicit = ellipse_ex2im(ellipse_transparent2ex(irisFitParams(ii,:)));
                fh=@(x,y) pFitImplicit(1).*x.^2 +pFitImplicit(2).*x.*y +pFitImplicit(3).*y.^2 +pFitImplicit(4).*x +pFitImplicit(5).*y +pFitImplicit(6);
                % superimpose the ellipse using fimplicit or ezplot (ezplot
                % is the fallback option for older Matlab versions)
                if exist('fimplicit','file')==2
                    fimplicit(fh,[1, videoSizeX, 1, videoSizeY],'Color', p.Results.irisColor,'LineWidth',1);
                    set(gca,'position',[0 0 1 1],'units','normalized')
                    axis off;
                else
                    plotHandle=ezplot(fh,[1, videoSizeX, 1, videoSizeY]);
                    set(plotHandle, 'Color', p.Results.irisColor)
                    set(plotHandle,'LineWidth',1);
                end
            end
        end
    end
    
    % add an instruction label
    if ~isempty(p.Results.controlFileName)
        instructionIdx = find ([instructions.frame] == ii);
        if ~isempty(instructionIdx)
            text_str = instructions(instructionIdx(end)).type;
            annotation('textbox',...
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
    
    % add the center of projection
    if ~isempty(p.Results.sceneGeometryFileName)
        plot(sceneGeometry.eyeCenter.X,sceneGeometry.eyeCenter.Y,['x' p.Results.sceneGeometryColor]);
    end
    
    % Save the frame to the scratch directory and close the figure
    tmp=getframe(frameFig);
    thisVideoFrame=uint8(tmp.cdata);
    thisFrameFileName= fullfile(scratchDirectoryName,sprintf('%05d.mat',ii));
    parsave(thisFrameFileName, thisVideoFrame);
    close(frameFig);
    
end

%% Save and cleanup

% Create a color map
cmap = [linspace(0,1,256)' linspace(0,1,256)' linspace(0,1,256)'];
cmap(1,:)=[1 0 0];
cmap(2,:)=[0 1 0];
cmap(3,:)=[0 0 1];
cmap(4,:)=[1 1 0];
cmap(5,:)=[0 1 1];
cmap(6,:)=[1 0 1];

% write the outputVideo to file
videoOutObj = VideoWriter(videoOutFileName,'Indexed AVI');
videoOutObj.FrameRate = p.Results.videoOutFrameRate;
videoOutObj.Colormap = cmap;
open(videoOutObj);

% loop through the frames and save them
for ii=1:nFrames
    thisFrameFileName=  fullfile(scratchDirectoryName,sprintf('%05d.mat',ii));
    tmp=load(thisFrameFileName);
    thisVideoFrame=tmp.dataout;
    indexedFrame = rgb2ind(thisVideoFrame, cmap, 'nodither');
    writeVideo(videoOutObj,indexedFrame);
end
% close the videoObj
clear videoOutObj

% delete the scratch directory
rmdir(scratchDirectoryName,'s');

% report completion of fit video generation
if strcmp(p.Results.verbosity,'full')
    toc
    fprintf('\n');
end

% Delete the parallel pool
if strcmp(p.Results.verbosity,'full')
    tic
    fprintf(['Closing parallel pool. Started ' char(datetime('now')) '\n']);
end
if p.Results.useParallel
    poolObj = gcp;
    if ~isempty(poolObj)
        delete(poolObj);
    end
end
if strcmp(p.Results.verbosity,'full')
    toc
    fprintf('\n');
end

end % function


%% LOCAL FUNCTIONS

function parsave(fname, dataout)
save(fname, 'dataout')
end
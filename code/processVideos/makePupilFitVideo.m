function makePupilFitVideo(videoInFileName, videoOutFileName, varargin)
% function makePupilFitVideo(perimeterVideoFileName, varargin)
%
% This routine will superimpose ellipses on frames of a video file.
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
%   glintColor
%   perimeterFileName
%   perimeterColor
%   ellipseFitFileName - full path to the .mat file that contains the
%       pupil tracking information.
%   ellipseColor
%   whichFieldToPlot - The name of the field within the ellipseFitData
%   	struct that is to be plotted
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
p.addParameter('tbtbRepoName','LiveTrackAnalysisToolbox',@ischar);

% Optional items to include in the video
p.addParameter('glintFileName',[],@(x)(isempty(x) | ischar(x)));
p.addParameter('glintColor','r',@ischar);
p.addParameter('perimeterFileName',[],@(x)(isempty(x) | ischar(x)));
p.addParameter('perimeterColor','w',@ischar);
p.addParameter('pupilFileName',[],@(x)(isempty(x) | ischar(x)));
p.addParameter('pupilColor','green',@ischar);
p.addParameter('whichFieldToPlot', [],@(x)(isempty(x) | ischar(x)));
p.addParameter('irisFileName',[],@(x)(isempty(x) | ischar(x)));
p.addParameter('irisCircleColor','magenta',@ischar);
p.addParameter('irisMaskColor','yellow',@ischar);
p.addParameter('controlFileName',[],@(x)(isempty(x) | ischar(x)));

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
else
    irisData=[];
end

% Read in and parse the control file if passed
if ~isempty(p.Results.controlFileName)
    instructions = importControlFile(p.Results.controlFileName);
else
    instructions(1).frame=[];
    instructions(1).type=[];
    instructions(1).params=[];
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
% read the video into memory, adjusting gamma if needed
for ii = 1:nFrames
    thisFrame = readFrame(videoInObj);
    sourceVideo(:,:,ii) = rgb2gray (thisFrame);
end
% close the video object
clear videoInObj

% prepare the outputVideo
outputVideo=zeros(videoSizeY,videoSizeX,3,nFrames,'uint8');


%% Loop through the frames
parfor (ii = 1:nFrames, nWorkers)

    % Update the progress display
    if strcmp(p.Results.verbosity,'full') && mod(ii,round(nFrames/50))==0
        fprintf('\b.\n');
    end
    
    % Create a figure
    frameFig = figure( 'Visible', 'off');

    % show the initial frame
    imshow(squeeze(sourceVideo(:,:,ii)), 'Border', 'tight');
    hold on
    
    % add glint
    if ~isempty(p.Results.glintFileName)
        if glintData.ellipseFittingError(ii)==1
            plot(glintData.X(ii),glintData.Y(ii),'*y');
        else
            plot(glintData.X(ii),glintData.Y(ii),['*' p.Results.glintColor]);
        end
    end
    
    % add pupil perimeter
    if ~isempty(p.Results.perimeterFileName)
        binP = squeeze(perimeter.data(:,:,ii));
        [Yp, Xp] = ind2sub(size(binP),find(binP));
        plot(Xp,Yp,['.' p.Results.perimeterColor], 'MarkerSize', 1);
    end
    
    % add irisMask
    if ~isempty(p.Results.irisFileName)
        binP = squeeze(irisData.mask(:,:,ii));
        boundaryPointsCell = bwboundaries(binP);
        if ~isempty(boundaryPointsCell)
            boundaryPoints=boundaryPointsCell{1};
            plot(boundaryPoints(:,2), boundaryPoints(:,1),['.' p.Results.irisMaskColor], 'MarkerSize', 1);
        end
    end
        
    % add ellipse fit
    if ~isempty(p.Results.pupilFileName)
        if sum(isnan(pupilFitParams(ii,:)))==0
            % build ellipse impicit equation
            pFitImplicit = ellipse_ex2im(ellipse_transparent2ex(pupilFitParams(ii,:)));
            fh=@(x,y) pFitImplicit(1).*x.^2 +pFitImplicit(2).*x.*y +pFitImplicit(3).*y.^2 +pFitImplicit(4).*x +pFitImplicit(5).*y +pFitImplicit(6);
            % superimpose the ellipse using fimplicit or ezplot
            if exist('fimplicit','file')==2
                fimplicit(fh,[1, videoSizeY, 1, videoSizeX],'Color', p.Results.pupilColor,'LineWidth',1.5);
                set(gca,'position',[0 0 1 1],'units','normalized')
                axis off;
            else
                plotHandle=ezplot(fh,[1, videoSizeY, 1, videoSizeX]);
                set(plotHandle, 'Color', p.Results.pupilColor)
            end
        end
    end
    
    % add iris circle fit
%     if ~isempty(p.Results.irisFileName)
%         if ~isnan(irisData.X(ii))
%             % build circle impicit equation
%             fh=@(x,y) (x-irisData.X(ii)).^2 +(y-irisData.Y(ii)).^2 - irisData.radius(ii).^2;
%             % superimpose the ellipse using fimplicit or ezplot
%             if exist('fimplicit','file')==2
%                 fimplicit(fh,[1, videoSizeY, 1, videoSizeX],'Color', p.Results.irisCircleColor,'LineWidth',1.5);
%                 set(gca,'position',[0 0 1 1],'units','normalized')
%                 axis off;
%             else
%                 plotHandle=ezplot(fh,[1, videoSizeY, 1, videoSizeX]);
%                 set(plotHandle, 'Color', p.Results.irisCircleColor)
%             end
%         end
%     end
    
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

    % Save the frame and close the figure
    tmp=getframe(frameFig);
    outputVideo(:,:,:,ii)=tmp.cdata;
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
    indexedFrame = rgb2ind(squeeze(outputVideo(:,:,:,ii)), cmap, 'nodither');
    writeVideo(videoOutObj,indexedFrame);
end
% close the videoObj
clear videoOutObj

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

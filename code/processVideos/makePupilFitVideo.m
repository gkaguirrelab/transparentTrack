function makePupilFitVideo(videoInFileName, videoOutFileName, varargin)
% function makePupilFitVideo(perimeterVideoFileName, varargin)
%
% This routine will superimpose ellipses on frames of a video file.
%
% INPUTS:
%   videoInFileName: full path to an .avi file. Points on the
%     boundary of the pupil should have a value of unity, and the frame
%     should be otherwise zero-filled. A frame that has no information
%     regarding the pupil (e.g., during a blink) should be zero-filled.
%   videoOutFileName: full path to an .avi file. Points on the
%     boundary of the pupil should have a value of unity, and the frame
%     should be otherwise zero-filled. A frame that has no information
%     regarding the pupil (e.g., during a blink) should be zero-filled.
%  ellipseFitDataFileName: full path to the .mat file that contains the
%     pupil tracking information.
%
% Optional key/value pairs (display and I/O)
%  'verbosity' - level of verbosity. [none, full]
%  'display' - controls a display of the fitting outcome. [none, full]
%  'videoOutFrameRate' - frame rate (in Hz) of saved video. Default 30.
%  'whichFieldToPlot' - The name of the field within the ellipseFitData
%     struct that is to be plotted
%
%  ellipseFitDataFileName

%% Parse vargin for options passed here
p = inputParser; p.KeepUnmatched = true;

% Required
p.addRequired('videoInFileName', @ischar);
p.addRequired('videoOutFileName', @ischar);

% Optional display and I/O params
p.addParameter('verbosity','none', @ischar);
p.addParameter('display','none', @ischar);
p.addParameter('videoOutFrameRate', 60, @isnumeric);
p.addParameter('nFrames',Inf,@isnumeric);

% Optional items to include in the video
p.addParameter('glintFileName',[],@(x)(isempty(x) | ischar(x)));
p.addParameter('glintColor','r',@isstring);
p.addParameter('perimeterFileName',[],@(x)(isempty(x) | ischar(x)));
p.addParameter('perimeterColor','w',@isstring);
p.addParameter('ellipseFitFileName',[],@(x)(isempty(x) | ischar(x)));
p.addParameter('ellipseColor','green',@isstring);
p.addParameter('whichFieldToPlot', [],@(x)(isempty(x) | ischar(x)));

% parse
p.parse(videoInFileName, videoOutFileName, varargin{:})

% Alert the user
if strcmp(p.Results.verbosity,'full')
    tic
    fprintf(['Creating and saving fit video. Started ' char(datetime('now')) '\n']);
    fprintf('| 0                      50                   100%% |\n');
    fprintf('.');
end

% Read in the glint file if passed
if ~isempty(p.Results.glintFileName)
    dataLoad = load(p.Results.glintFileName);
    glintData = dataLoad.glintData;
    clear dataLoad
end

% Read in the perimeter file if passed
if ~isempty(p.Results.perimeterFileName)
    dataLoad = load(p.Results.perimeterFileName);
    perimeter = dataLoad.perimeter;
    clear dataLoad
end

% Read in the ellipseFitData file if passed
if ~isempty(p.Results.ellipseFitFileName)
    dataLoad = load(p.Results.ellipseFitFileName);
    ellipseFitData = dataLoad.ellipseFitData;    
    clear dataLoad
    ellipseFitParams = ellipseFitData.(p.Results.whichFieldToPlot);
end

% Create a figure
if strcmp(p.Results.display,'full')
    frameFig = figure( 'Visible', 'on');
else
    frameFig = figure( 'Visible', 'off');
end


% open inVideoObject
inVideoObj = VideoReader(videoInFileName);
% get number of frames
if p.Results.nFrames == Inf
    nFrames = floor(inVideoObj.Duration*inVideoObj.FrameRate);
else
    nFrames = p.Results.nFrames;
end

% get video dimensions
videoSizeX = inVideoObj.Width;
videoSizeY = inVideoObj.Height;

% Open the video out object
outVideoObj = VideoWriter(videoOutFileName);
outVideoObj.FrameRate = p.Results.videoOutFrameRate;
open(outVideoObj);


% Loop through the frames
for ii=1:nFrames
    
    % Update the progress display
    if strcmp(p.Results.verbosity,'full') && mod(ii,round(nFrames/50))==0
            fprintf('.');
    end
    
    % Read the frame, adjust gamma, make gray
    thisFrame = readFrame(inVideoObj);
    thisFrame = rgb2gray (thisFrame);

    % show the initial frame
        imshow(thisFrame, 'Border', 'tight');
    hold on
    
    % add glint
    if ~isempty(p.Results.glintFileName)
        plot(glintData.X(ii),glintData.Y(ii),['*' p.Results.glintColor]);
    end
    
    % add pupil perimeter
    if ~isempty(p.Results.perimeterFileName)
        binP = squeeze(perimeter.data(:,:,ii));
        [Yp, Xp] = ind2sub(size(binP),find(binP));
        plot(Xp,Yp,['.' p.Results.perimeterColor], 'MarkerSize', 1);
    end
    
    % add ellipse fit
    if ~isempty(p.Results.ellipseFitFileName)
        if sum(isnan(ellipseFitParams(ii,:)))==0
            % build ellipse impicit equation
            pFitImplicit = ellipse_ex2im(ellipse_transparent2ex(ellipseFitParams(ii,:)));
            fh=@(x,y) pFitImplicit(1).*x.^2 +pFitImplicit(2).*x.*y +pFitImplicit(3).*y.^2 +pFitImplicit(4).*x +pFitImplicit(5).*y +pFitImplicit(6);
            
            % superimpose the ellipse using fimplicit
            hold on
            if strcmp(version('-release'),'2016a')
                ezplot(fh,[1, videoSizeY, 1, videoSizeX]);
            else
                fimplicit(fh,[1, videoSizeY, 1, videoSizeX],'Color', p.Results.ellipseColor);
                set(gca,'position',[0 0 1 1],'units','normalized')
                axis off;
            end
        end
    end
    
    % Write the frame to the file
    writeVideo(outVideoObj,getframe(frameFig));
    
    hold off
end

% close the video object
close(outVideoObj);

% close the figure
close(frameFig)

% report completion of fit video generation
if strcmp(p.Results.verbosity,'full')
    fprintf('\n');
    toc
end


end


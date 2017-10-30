function resizeAndCropVideo(inputVideoName, outputVideoName, varargin)
% function resizeAndCropVideo(inputVideoName, outputVideoName, varargin)
%
%  This fuction crops and resizes a video according to the specified
%  parameters. Using the default options the routine will scale and crop a
%  VGA video to match the LiveTrack standard size video. The
%  video is also converted to gray, if requested.

% Output
%   an AVI video is saved out.
%
% Input (required)
%	inputVideoName - full path to the video to deinterlace
%	outputVideoName - full path to the deinterlaced output video
%
% Options (analysis)
%   resizeVideo - [Y X] desired output video resolution. (keep default to
%       get livetrack format)
%   cropVideo - [firstX firstY lastX lastY] position of first and last
%           pixels to include in the crop. (keep default to get livetrack
%           format)
%   keepOriginalSize - option to skip video resizing.
%
% Options (verbosity and display)
%   verbosity - controls console status updates
%
% Options (flow control)
%  nFrames' - analyze fewer than the total number of frames.
%  startFrame - which frame to start on
%
%
% Usage examples
%  resizeAndCropVideo(inputVideoName,outputVideoName);
%  resizeAndCropVideo(inputVideoName,outputVideoName, 'nFrames', 1000) % this will
%       process just the first 1000 frames of the video


%% parse input and define variables

p = inputParser;
% required input
p.addRequired('inputVideoName',@isstr);
p.addRequired('outputVideoName',@isstr);

% optional inputs
p.addParameter('resizeVideo',[486 720]/2, @isnumeric);
p.addParameter('cropVideo', [1 1 319 239], @isnumeric);
p.addParameter('keepOriginalSize', false, @islogic);
p.addParameter('convertToGray',true,@islogical)

% verbosity
p.addParameter('verbosity', 'none', @isstr);

% flow control
p.addParameter('nFrames',Inf,@isnumeric);
p.addParameter('startFrame',1,@isnumeric);


%parse
p.parse(inputVideoName,outputVideoName,varargin{:})

% define variables
resizeVideo = p.Results.resizeVideo;
cropVideo = p.Results.cropVideo;
keepOriginalSize = p.Results.keepOriginalSize;

%% Prepare Video

% load video
inObj = VideoReader(inputVideoName);

% create outputVideo object
outObj = VideoWriter(outputVideoName);
outObj.FrameRate = inObj.FrameRate;
open(outObj);

% option to manually set numFrames
if p.Results.nFrames == Inf
    nFrames = floor(inObj.Duration*inObj.FrameRate);
else
    nFrames=p.Results.nFrames;
end

% alert the user
if strcmp(p.Results.verbosity,'full')
    tic
    fprintf(['Resizing and cropping video. Started ' char(datetime('now')) '\n']);
    fprintf('| 0                      50                   100%% |\n');
    fprintf('.');
end

% Resize and crop, save
for ii = p.Results.startFrame:nFrames
    % increment the progress bar
    if strcmp(p.Results.verbosity,'full') && mod(ii,round(nFrames/50))==0
        fprintf('.');
    end
    thisFrame = readFrame(inObj);
    if p.Results.convertToGray
        thisFrame = rgb2gray(thisFrame);
    end
    if keepOriginalSize == 0
        tmp = imresize(thisFrame,resizeVideo);
        thisFrame = imcrop(tmp,cropVideo);
    end
    writeVideo(outObj,thisFrame);
end

clear inObj outObj

% report completion of analysis
if strcmp(p.Results.verbosity,'full')
    fprintf('\n');
    toc
    fprintf('\n');
end

end % function
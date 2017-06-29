function resizeAndCropVideo(inputVideoName, outputVideoName, varargin)
% function resizeAndCropVideo(inputVideoName, outputVideoName, varargin)
%
%  This fuction converts the video to a "gray frames array" that is stored
%  in the memory and ready to be tracked or written to file. With the
%  default options the routine will scale and crop the video to livetrack
%  standard size.

% Output
% ======
%       a gray video is saved out
%
% Input
% =====
%       inputVideoName
%       outputVideoName
%
% Options
% =======
%       numberOfFrames : number of frames to process. If not specified or
%           Inf will process the full video.
%       resizeVideo : [Y X] desired output video resolution. (recommended: keep default)
%       cropVideo : [firstX firstY lastX lastY] position of first and last
%           pixels to include in the crop. (recommended: keep default)
%       keepOriginalSize : option to skip video resizing.
%
%  NOTE: if processing videos acquired with the LiveTrack+V.TOP hardware
%  setup, do not alter the default resizing and cropping video options
%
%
% Usage examples
% ==============
%  resizeAndCropVideo(inputVideoName,outputVideoName);
%  resizeAndCropVideo(inputVideoName,outputVideoName, 'numberOfFrames', 1000) % this will
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
p.addParameter('nFrames', Inf, @isnumeric);
p.addParameter('verbosity', 'none', @ischar);


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

if strcmp(p.Results.verbosity,'full')
    tic
    fprintf(['Resizing and cropping video. Started ' char(datetime('now')) '\n']);
    fprintf('| 0                      50                   100%% |\n');
    fprintf('.');
end

% Resize and crop, save
for ii = 1:nFrames
    if strcmp(p.Results.verbosity,'full') && mod(ii,round(nFrames/50))==0
        fprintf('.');
    end
    thisFrame = readFrame(inObj);
    tmp = rgb2gray(thisFrame);
    if keepOriginalSize == 0
        tmp2 = imresize(tmp,resizeVideo);
        tmp = imcrop(tmp2,cropVideo);
    end
    writeVideo(outObj,tmp);
    % increment progress bar
end

clear inObj outObj

% report completion of analysis
if strcmp(p.Results.verbosity,'full')
    fprintf('\n');
    toc
    fprintf('\n');
end

end % function
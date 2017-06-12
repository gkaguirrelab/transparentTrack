function [grayI] = prepareVideo(inputVideo, varargin)

%  This fuction converts the video to a "gray frames array" that is stored
%  in the memory and ready to be tracked or written to file. With the
%  default options the routine will scale and crop the video to livetrack
%  standard size.

% Output
% ======
%       grayI = 3D array of cropped and grey scaled video frames.
% 
% Input
% =====
%       inputVideo
% 
% Options
% =======
%       numberOfFrames : number of frames to process. If not specified will
%           process the full video.
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
%  [grayI] = prepareVideo(inputVideo);
%  [grayI] = prepareVideo(inputVideo, 'numberOfFrames', 1000) % this will
%       process just the first 1000 frames of the video


%% parse input and define variables

p = inputParser;
% required input
p.addRequired('inputVideo',@isstr);
% optional inputs
defaultResize = [486 720]/2;
defaultCrop = [1 1 319 239];
defaultNumFrames = inf;
defaultKeepOriginalSize = false;
p.addParameter('resizeVideo', defaultResize, @isnumeric);
p.addParameter('cropVideo', defaultCrop, @isnumeric);
p.addParameter('numberOfFrames', defaultNumFrames, @isnumeric);
p.addParameter('keepOriginalSize', defaultKeepOriginalSize, @islogic);
%parse
p.parse(inputVideo,varargin{:})

% define variables
resizeVideo = p.Results.resizeVideo;
cropVideo = p.Results.cropVideo;
numberOfFrames = p.Results.numberOfFrames;
keepOriginalSize = p.Results.keepOriginalSize;

%% Prepare Video

% load video
disp('Loading video file...');
inObj = VideoReader(inputVideo);

% option to manually set numFrames
if numberOfFrames ~= Inf
    numFrames = numberOfFrames;
else
    numFrames = floor(inObj.Duration*inObj.FrameRate);
end

% initialize gray image array for default format (for faster processing)
if resizeVideo == defaultResize
    grayI = zeros([240 320 numFrames],'uint8');
else
    warning('Using non default resizing option. The gray image array will not be inizialized (converting will take longer)')
end

% Convert to gray, resize, crop
disp('Converting video to gray frames array, may take a while...');
for i = 1:numFrames
    thisFrame           = readFrame(inObj);
    tmp                 = rgb2gray(thisFrame);
    if keepOriginalSize == 0
        tmp2 = imresize(tmp,resizeVideo);
        tmp = imcrop(tmp2,cropVideo);
    end
    grayI(:,:,i) = tmp;
end

clear RGB inObj

disp('> done.');
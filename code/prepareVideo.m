function [grayI, params] = prepareVideo(params)

%  This fuction converts the video to a "gray frames array" that are store
%  in the memory and ready to be tracked or written to file. As default,
%  the routine will crop the video to livetrack size.

% Output : 
%       grayI = 4D array of cropped and grey scale video frames.
% 
% Input:
%       params.inVideo
%       params.keepOriginalSize
%       params.imageSize
%       params.imageCrop
%       params.forceNumFrames


%% set default params
if ~isfield (params, 'keepOriginalSize')
    params.keepOriginalSize = 0;
end
if ~isfield(params,'imageSize')
    params.imageSize = [486 720]/2;
end
if ~isfield(params,'imageCrop')
    params.imageCrop = [1 1 319 239];
end

%% Load video
disp('Loading video file...');
inObj = VideoReader(params.inVideo);
numFrames = floor(inObj.Duration*inObj.FrameRate);

% option to overwrite numFrames (for processing small sections of video)
if isfield(params,'forceNumFrames')
    numFrames = params.forceNumFrames;
end

% initialize gray image array
grayI                   = zeros([240 320 numFrames],'uint8');

disp('Converting video to gray frames array, may take a while...');
% Convert to gray, resize, crop to livetrack size
for i = 1:numFrames
    thisFrame           = readFrame(inObj);
    tmp                 = rgb2gray(thisFrame);
    if params.keepOriginalSize == 0
        tmp2 = imresize(tmp,params.imageSize);
        tmp = imcrop(tmp2,params.imageCrop);
    end
    grayI(:,:,i) = tmp;
end

clear RGB inObj

disp('> done.');
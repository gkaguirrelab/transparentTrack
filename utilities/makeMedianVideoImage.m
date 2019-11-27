function [imageMedian, imageSD] = makeMedianVideoImage(videoInFileName, varargin)
% Creates the median video image, masking the pupil and iris
%
% Syntax:
%  makeMedianVideoImage(videoInFileName)
%
% Description:
%   foo
%
% Inputs:
%	videoInFileName       - Full path to an .avi file; typically the "gray"
%
% Optional key/value pairs (display and I/O):
%
% Outputs:
%   None
%

%% Parse vargin for options passed here
p = inputParser; p.KeepUnmatched = false;

% Required
p.addRequired('videoInFileName', @ischar);

% Optional flow control params
p.addParameter('startFrame', 1, @isnumeric);
p.addParameter('nFrames', Inf, @isnumeric);

% parse
p.parse(videoInFileName, varargin{:})


% Prepare the video object
videoInObj = videoIOWrapper(videoInFileName,'ioAction','read');

% get number of frames
if p.Results.nFrames == Inf
    nFrames = floor(videoInObj.Duration*videoInObj.FrameRate);
else
    nFrames = p.Results.nFrames;
end

% get start and end times
endTime = (p.Results.startFrame+nFrames) / videoInObj.FrameRate;
startTime = p.Results.startFrame / videoInObj.FrameRate;

% get video dimensions
videoSizeX = videoInObj.Width;
videoSizeY = videoInObj.Height;

sourceFrames = zeros(1,videoSizeY,videoSizeX);

% Make sure that the start and end times are in the right order
if startTime > endTime
    tmp = endTime;
    endTime = startTime;
    startTime = tmp;    
end

%% Loop through the frames
% Obtain the median for each second of video
for ii = startTime:(1/videoInObj.FrameRate):endTime
    % read the source video frame into memory
    videoInObj.CurrentTime=ii;
    frame = rgb2gray(readFrame(videoInObj));

    sourceFrames(end+1,:,:)=frame;
end

% Derive the median and SD images
sourceFrames = sourceFrames(2:end,:,:);
imageMedian = squeeze(median(sourceFrames,1));
imageSD = squeeze(std(sourceFrames,1));

% close the video object
clear videoInObj


end % function

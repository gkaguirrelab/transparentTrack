function [imageMedian, imageSD] = makeMedianVideoImage(videoInFileName, varargin)
% Returns the median image for a specified set of video frames by index
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
%   startFrame            - The first frame of the median, indexed from 1
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
startFrame = p.Results.startFrame;
endFrame = startFrame+nFrames-1;

% get video dimensions
videoSizeX = videoInObj.Width;
videoSizeY = videoInObj.Height;

% A variable to hold the frames as we read them
sourceFrames = zeros(1,videoSizeY,videoSizeX);

% Read up to the startFrame
curFrame = 1;
while curFrame < startFrame
    readFrame(videoInObj);
    curFrame = curFrame + 1;
end

%% Loop through the frames
for ii = startFrame:endFrame
    
    % read the source video frame into memory
    frame = rgb2gray(readFrame(videoInObj));
    
    % store this frame
    sourceFrames(end+1,:,:)=frame;
end

% Derive the median and SD images
imageMedian = squeeze(median(sourceFrames,1));
imageSD = squeeze(std(sourceFrames,1));

% close the video object
clear videoInObj


end % function

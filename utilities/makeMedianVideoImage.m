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
p = inputParser; p.KeepUnmatched = true;

% Required
p.addRequired('videoInFileName', @ischar);

% Optional flow control params
p.addParameter('startFrame', 1, @isnumeric);
p.addParameter('nFrames', Inf, @isnumeric);

% Optional video items
p.addParameter('chunkSizeSecs', 20, @isscalar);
p.addParameter('outputFileName', [], @(x)(isempty(x) | ischar(x)));

% parse
p.parse(videoInFileName, varargin{:})

chunkSizeSecs = p.Results.chunkSizeSecs;


% Prepare the video object
videoInObj = videoIOWrapper(videoInFileName,'ioAction','read');

% get number of frames
if p.Results.nFrames == Inf
    nFrames = floor(videoInObj.Duration*videoInObj.FrameRate);
else
    nFrames = p.Results.nFrames;
end

% get start and end times
startTime = (ceil(p.Results.startFrame / videoInObj.FrameRate)/chunkSizeSecs)*chunkSizeSecs;
endTime = floor(((p.Results.startFrame+nFrames) / videoInObj.FrameRate) /chunkSizeSecs)*chunkSizeSecs;

% get video dimensions
videoSizeX = videoInObj.Width;
videoSizeY = videoInObj.Height;

sourceFrames = zeros(1,videoSizeY,videoSizeX);

%% Loop through the frames
% Obtain the median for each second of video
for ii = startTime:chunkSizeSecs:endTime
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

% Save the images if requested
if ~isempty(p.Results.outputFileName)
    save(p.Results.outputFileName,'imageMedian','imageSD');
end


end % function

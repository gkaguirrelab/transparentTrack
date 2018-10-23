function [imageMedian, imageSD] = makeMedianVideoImage(videoInFileName, varargin)
% Create and store a video that displays the results of eye tracking
%
% Syntax:
%  makeMedianVideoImage(videoInFileName)
%
% Description:
%   This routine creates an integrated fit video that illustrates the
%   position of the pupil and glint and indicates any control instructions.
%
% Inputs:
%	videoInFileName       - Full path to an .avi file; typically the "gray"
%
% Optional key/value pairs (display and I/O):
%  'pupilFileName'        - Full path to a pupil file used to remove the
%                           pupil and iris from the image.
%  'fitLabel'             - The field of the pupilData file that contains
%
% Outputs:
%   None
%

%% Parse vargin for options passed here
p = inputParser; p.KeepUnmatched = true;

% Required
p.addRequired('videoInFileName', @ischar);

% Optional flow control params
p.addParameter('nFrames', Inf, @isnumeric);

% Optional video items
p.addParameter('pupilFileName', [], @(x)(isempty(x) | ischar(x)));
p.addParameter('fitLabel', 'initial', @(x)(isempty(x) | ischar(x)));
p.addParameter('chunkSizeSecs', 5, @isscalar);

% parse
p.parse(videoInFileName, varargin{:})

chunkSizeSecs = p.Results.chunkSizeSecs;

% Read in the pupilData file if passed
if ~isempty(p.Results.pupilFileName)
    dataLoad = load(p.Results.pupilFileName);
    pupilData = dataLoad.pupilData;
    clear dataLoad
    ellipseFitParams = pupilData.(p.Results.fitLabel).ellipses.values;
    ellipseFitRMSE = pupilData.(p.Results.fitLabel).ellipses.RMSE;
    if isfield(pupilData.(p.Results.fitLabel),'eyePoses')
        eyePoses = pupilData.(p.Results.fitLabel).eyePoses.values;
    else
        eyePoses = [];
    end
else
    ellipseFitParams=[];
    ellipseFitRMSE=[];
    eyePoses=[];
end

% Prepare the video object
videoInObj = videoIOWrapper(videoInFileName,'ioAction','read');

% get number of frames
if p.Results.nFrames == Inf
    nFrames = floor(videoInObj.Duration*videoInObj.FrameRate);
else
    nFrames = p.Results.nFrames;
end

% get video dimensions
videoSizeX = videoInObj.Width;
videoSizeY = videoInObj.Height;

sourceFrames = zeros(1,videoSizeY,videoSizeX);

%% Loop through the frames
% Obtain the median for each second of video
for ii = chunkSizeSecs:chunkSizeSecs:floor(videoInObj.Duration/chunkSizeSecs)*chunkSizeSecs
    % read the source video frame into memory
    videoInObj.CurrentTime=ii;
    frame = rgb2gray(readFrame(videoInObj));

    % If the pupil ellipse is available, nan out an area around the center
    % of the pupil
    if ~isempty(p.Results.pupilFileName)
        if ~isempty(ellipseFitParams)
            if sum(isnan(ellipseFitParams(ii,:)))==0
                frame = double(rgb2gray(insertShape(frame,'FilledCircle',[ellipseFitParams(ii,1) ellipseFitParams(ii,2) 100],'Color','black','Opacity',1)));
            end
        end
    end    
    frame=double(frame);
    frame(frame==0)=nan;

    sourceFrames(end+1,:,:)=frame;
end

% Derive the median and SD images
sourceFrames = sourceFrames(2:end,:,:);
imageMedian = squeeze(median(sourceFrames,1));
imageSD = squeeze(std(sourceFrames,1));

% close the video object
clear videoInObj


end % function

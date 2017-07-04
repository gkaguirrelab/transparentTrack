function blinkFrames = findBlinks(glintFile, varargin)

% this function flags the blinks based on the glint location and pupil
% perimeter extraction. User can specify the amount of frames around the
% actual blink to be flagged as a blink as well (i.e. the blink window).
% 
% Input params
% ============
%       glintFile : matfile with glint tracking data
%       perimeterVideo : the perimeter video.
%       
% Options
% =======
%        extendBlinkWindow : number of additional frames flagged as blink
%           before and after the detected blink frame. If more than one
%           consecutive blink frame is detected, the extra frames will be
%           added BEFORE and AFTER the entire blink sequence. (default = 0)
% 
% Output
% ======
%       blinkFrames : index of the video frames flagged as blinks.
% 
% 
% Usage examples
% ==============
% 
%  blinkWindow = 5;
%  findBlinks(glintFile, perimeterVideo, controlFile, 'extendBlinkWindow', extendBlinkWindow)
% 
%% parse input and define variables

p = inputParser;
% required input
p.addRequired('glintFile',@isstr);

% optional inputs
extendBlinkWindowDefault = 0;
p.addParameter('extendBlinkWindow', extendBlinkWindowDefault, @isnumeric);

%parse
p.parse(glintFile, varargin{:})

% define optional variables values
extendBlinkWindow = p.Results.extendBlinkWindow;

%% load data

% glint
load(glintFile)

%% find the blinks

% since the glint is a reflection of a fixed light source, its location
% will be fairly constant thoroughout the video. During blinks, the routine
% may pick up a bright dot on the eyelight and wrongly identfy it as a
% legitimate glint. This step will remove all glint tracked in
% non-plausible location, when compared to the average glint location in
% the video.

% locate all nans
blinks = find (isnan(glintData.X));

% locate all points further away than 2 std from the glint mean
farX = find(glintData.X > nanmean(glintData.X) + (2 * nanstd(glintData.X)) | glintData.X < nanmean(glintData.X) - (2 * nanstd(glintData.X)));
farY = find(glintData.Y > nanmean(glintData.Y) + (2 * nanstd(glintData.Y)) | glintData.Y < nanmean(glintData.Y) - (2* nanstd(glintData.Y)));

farGlints = union(farX, farY);

% combine to obtain all blinks
blinkFrames = union(blinks, farGlints);


% %% extend blink window === TO DEVELOP
% 
% % find blink windows
% blinkWindows = find (diff(blinks) > 1);
% 
% for bb = 1:length(blinks)
%     
%

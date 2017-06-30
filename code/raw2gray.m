function raw2gray(inputVideoName,outputVideoName,varargin)
% function raw2gray(inputVideoName,outputVideoName,varargin)
%
% This function converts the raw, RBG 30Hz video acquired with the V.top
% device into a 60Hz, grayscale video. Additionally, the video is cropped
% to correspond to the same field of view as that returned by the LiveTrack
% device.
% 
% This process happens in two steps:
% 1. the original video is deinterlaced using the function deinterlaceVideo
% 2. the resulting video is converted to gray and resized to LiveTrack
% format using the function resizeAndCropVideo.
% 
% Note that between step 1 and 2 a large temporary video file is created,
% but it is deleted at the end of the routine, as it won't be necessary in
% the subsequent eyetracking steps.
% 
% 
%% parse input
p = inputParser;

% required input
p.addRequired('inputVideoName',@isstr);
p.addRequired('outputVideoName',@isstr);

% optional inputs
p.addParameter('nFrames', Inf, @isnumeric);
p.addParameter('verbosity', 'none', @ischar);

% parse
p.parse(inputVideoName,outputVideoName,varargin{:});


%% deinterlace video

% get file parts for inputVideoName
[~,nameRoot,~] = fileparts(inputVideoName);

% Check if the name root ends in a '_raw' suffix (standard saving format
% for V.TOP acquisitons). If so remove suffix to get the real nameRoot.
if regexp(nameRoot, regexptranslate('wildcard','*_raw'))
    nameRoot = nameRoot(1:end-4); %runs
end

% get saving path from outputVideoName
[pathstr,~,~] = fileparts(outputVideoName);

% set deinterlaced video name (large temporary file)
deinterlacedVideoName = fullfile(pathstr,[nameRoot '_60hz.avi']);

% check if deinterlaced video already exists
if ~exist(deinterlacedVideoName, 'file')
    deinterlaceVideo(inputVideoName, deinterlacedVideoName, 'nFrames', p.Results.nFrames, 'verbosity', p.Results.verbosity);
end

% adjust the video to LiveTrack format
resizeAndCropVideo(deinterlacedVideoName, outputVideoName, 'nFrames', p.Results.nFrames, 'verbosity', p.Results.verbosity);

%% remove the temporary 60Hz video
delete (deinterlacedVideoName)


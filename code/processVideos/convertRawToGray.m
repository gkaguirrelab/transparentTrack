function convertRawToGray(inputVideoName,outputVideoName,varargin)
% convertRawToGray(inputVideoName,outputVideoName,varargin)
%
% This routine is idiosyncratic to the analysis of interlaced, analog
% videos. It converts raw, RBG 30Hz video  into a 60Hz, grayscale video.
% Additionally, the video is cropped. The implementation here is fairly
% specific to the GKAguirreLab at U Penn and data collected using the
% Cambridge Research Systems LiveTrack device.
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
% This code could be modified for use on other videos to skip or modify one
% or both components

%% parse input
p = inputParser; p.KeepUnmatched = true;

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

end % main function



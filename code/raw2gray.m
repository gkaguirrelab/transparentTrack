function raw2gray(inputVideoName,outputVideoName,varargin)

% raw2gray turns the raw RBG 30Hz video acquired with the V.top device in a
% 60Hz gray color video in the same format as the one acquired by the
% LiveTrack device.
% 
% This process happens in two steps:
% 1. the original video is deinterlaced using the function deinterlaceVideo
% 2. the resulting video is converted to gray and resized to LiveTrack
% format using the function prepareVideo.
% 
% Note that between step 1 and 2 a large temporary video file is create,
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
p.addParameter('numberOfFrames', Inf, @isnumeric);

% parse
p.parse(inputVideoName,outputVideoName,varargin);

% define variables
numberOfFrames = p.Results.numberOfFrames;

%% deinterlace video

% get file parts for inputVideoName
[pathstr,name,~] = fileparts(inputVideoName);

% set deinterlaced video name (large temporary file)
deinterlacedVideoName = fullfile(pathstr,[name '_60Hz.avi']);

% check if deinterlaced video already exists
if ~exist(deinterlacedVideoName, 'file')
    deinterlaceVideo(inputVideoName, deinterlacedVideoName);
end

%% get the video in livetrack format
prepareVideo(deinterlacedVideoName, outputVideoName,'numberOfFrames',numberOfFrames);

%% remove the temporary 60Hz video
delete (deinterlacedVideoName)


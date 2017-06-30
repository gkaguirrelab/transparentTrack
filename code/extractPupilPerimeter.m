function [perimeterParams] = extractPupilPerimeter(grayVideoName, perimeterVideoName,varargin)
% 
% This function thresholds the video to extract the pupil perimeter and
% saves out a BW video showing the pupil perimeter only.
% 
% 
% Input params
% ============
%       grayVideoName : name and path of the gray video to track
%       perimeterVideoName : name of the output avi in which to save the output.
%       
% Options
% =======
%       frameRate: frame rate for the output video (default 60)
%       gammaCorrection : gamma correction to be applied in current frame 
%           (default 1, typical range [0.5 1.8])
%       pupilCircleThresh : threshold value to locate the glint for circle
%           fitting (default 0.6, typical range [0.5 0.9])
%       pupilRange : initial radius range for circle fitting of the glint
%       (default [30 90]). This value gets dynamically updated.
%       pupilEllipseThresh : threshold value to locate the glint for
%           ellipse fitting (default 0.9, typical range [0.8 0.98])
% 
% Output
% ======
%       perimeterVideo, perimeterParams.
% 
% 
% Usage examples
% ==============
% 
%  gammaCorrection = 1.2;
%  pupilEllipseThresh = 0.93;
%  perimeterParams = extractPupilPerimeter(grayVideoName, perimeterVideoName, 'gammaCorrection', gammaCorrection, 'pupilEllipseThresh', pupilEllipseThresh)
% 
%% parse input and define variables

p = inputParser;
% required input
p.addRequired('grayVideoName',@isstr);
p.addRequired('perimeterVideo',@isstr);

% optional inputs
frameRateDefault = 60;
gammaCorrectionDefault = 1;
pupilCircleThreshDefault =  0.06;
pupilRangeDefault = [30 90];
pupilEllipseThreshDefault = 0.95;
p.addParameter('frameRate', frameRateDefault, @isnumeric);
p.addParameter('gammaCorrection', gammaCorrectionDefault, @isnumeric);
p.addParameter('pupilCircleThresh', pupilCircleThreshDefault, @isnumeric);
p.addParameter('pupilRange', pupilRangeDefault, @isnumeric);
p.addParameter('pupilEllipseThresh', pupilEllipseThreshDefault, @isnumeric);

% Optional display and I/O params
p.addParameter('verbosity','none',@ischar);

% Environment parameters
p.addParameter('tbSnapshot',[],@(x)(isempty(x) | isstruct(x)));
p.addParameter('timestamp',char(datetime('now')),@ischar);
p.addParameter('username',char(java.lang.System.getProperty('user.name')),@ischar);
p.addParameter('hostname',char(java.net.InetAddress.getLocalHost.getHostName),@ischar);

% parse
p.parse(grayVideoName, perimeterVideoName, varargin{:})

% define optional variables values
frameRate = p.Results.frameRate;
gammaCorrection = p.Results.gammaCorrection;
pupilCircleThresh =  p.Results.pupilCircleThresh;
pupilRange = p.Results.pupilRange;
pupilEllipseThresh = p.Results.pupilEllipseThresh;
verbosity =  p.Results.verbosity;



%% read input video file
% load pupilPerimeter
inObj = VideoReader(grayVideoName);

% get number of frames
numFrames = floor(inObj.Duration*inObj.FrameRate);

%% initiate output video object

outObj = VideoWriter(perimeterVideoName);
outObj.FrameRate = frameRate;
open(outObj);

%% extract pupil perimeter

% initiate progress bar 
progBar = ProgressBar(numFrames,'Extracting perimeter...');

% open a figure
ih = figure;

% loop through gray frames
for ii = 1:numFrames
    
    % Get the frame
    I = readFrame(inObj);
    
     % adjust gamma for this frame
    I = imadjust(I,[],[],gammaCorrection);
    
    I = rgb2gray (I);
    
    % Show the frame
    imshow(I, 'Border', 'tight');

    % track with circles 
    % feed default values for glint (glint data won't be stored now)
    glintCircleThresh = 0.999;
    glintRange = [10 30];
    
    [pCenters, pRadii,pMetric,~,~,~, pupilRange, ~] = circleFit(I,pupilCircleThresh,glintCircleThresh,pupilRange,glintRange);
    
    if isempty(pCenters) %no pupil circle patch was found
        % make the frame black and save it
        I = zeros(size(I));
        imshow(I, 'Border', 'tight');
        frame   = getframe(ih);
        writeVideo(outObj,frame);
        % increment progress bar
        if ~mod(ii,10);progBar(ii);end;
        continue
    else
        % get pupil perimeter
        [binP] = getPupilPerimeter(I,pCenters,pRadii, pupilEllipseThresh);
        
        % convert BW frame
        I = im2uint8(binP);
        
        % Show the frame
        imshow(I, 'Border', 'tight');


        % write frame to thresholded video
        frame   = getframe(ih);
        writeVideo(outObj,frame);
        
        % increment progress bar
        if ~mod(ii,10);progBar(ii);end;
    end
end % loop through gray frames

%% close video
close(ih);
clear outObj

%% export perimeter params
perimeterParams = p.Results;

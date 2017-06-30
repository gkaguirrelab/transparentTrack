function [perimeter] = extractPupilPerimeter(grayVideoName, perimeterVideoName,varargin)
% function [perimeterParams] = extractPupilPerimeter(grayVideoName, perimeterVideoName,varargin)
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
%           fitting (default 0.06, typical range [0.04 0.09])
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
p.addParameter('glintCircleThresh', 0.999, @isnumeric);
p.addParameter('glintRange', [10 30], @isnumeric);

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
glintCircleThresh = p.Results.glintCircleThresh;
glintRange = p.Results.glintCircleThresh;

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

% alert the user
if strcmp(p.Results.verbosity,'full')
    tic
    fprintf(['Tracking the glint. Started ' char(datetime('now')) '\n']);
    fprintf('| 0                      50                   100%% |\n');
    fprintf('.');
end

% open a figure
ih = figure;

% loop through gray frames
for ii = 1:numFrames
    % increment the progress bar
    if strcmp(p.Results.verbosity,'full') && mod(ii,round(nFrames/50))==0
        fprintf('.');
    end
    
    % Get the frame
    I = readFrame(inObj);
    
     % adjust gamma for this frame
    I = imadjust(I,[],[],gammaCorrection);
    
    I = rgb2gray (I);
    
    % Show the frame
    imshow(I, 'Border', 'tight');

    % track with circles 
    [pCenters, pRadii,~,~,~,~, pupilRange, ~] = circleFit(I,pupilCircleThresh,glintCircleThresh,pupilRange,glintRange);
    
    if isempty(pCenters) %no pupil circle patch was found
        % make the frame black and save it
        I = zeros(size(I));
        imshow(I, 'Border', 'tight');
        frame   = getframe(ih);
        writeVideo(outObj,frame);
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
        
    end
end % loop through gray frames

%% close video
close(ih);
clear outObj

%% save mat file with analysis details and save it
perimeter.meta = p.Results;

[perimeterMatPath, perimeterMatName, ~] = fileparts(perimeterVideoName);

save (fullfile(perimeterMatPath,[perimeterMatName '.mat']),'perimeter')

% report completion of analysis
if strcmp(p.Results.verbosity,'full')
    fprintf('\n');
    toc
    fprintf('\n');
end

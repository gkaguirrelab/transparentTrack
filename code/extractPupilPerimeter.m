function perimeterParams = extractPupilPerimeter(grayI, perimeterVideoPath,varargin)

% This function thresholds the video to extract the pupil perimeter and
% saves out a BW video showing the pupil perimeter only.
% 

% Input params
% ============
%       grayI : 3D array of gray frames to track
%       perimeterVideoPath : name of the output avi in which to save the output.
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
%  define options beforehand:

%  gammaCorrection = 1.2;
%  pupilEllipseThresh = 0.93;
%  perimeterParams = extractPupilPerimeter(grayI, perimeterVideoPath, 'gammaCorrection', gammaCorrection, 'pupilEllipseThresh', pupilEllipseThresh)
% 
%  OR
% 
%  define options directly:
% 
%  perimeterParams = extractPupilPerimeter(grayI, perimeterVideoPath, 'gammaCorrection', 1.2, 'pupilEllipseThresh', 0.93)


%% parse input and define variables

p = inputParser;
% required input
p.addRequired('grayI');
p.addRequired('perimeterVideoPath',@isstr);

% optional inputs
frameRateDefault = 60;
gammaCorrectionDefault = 1;
pupilCircleThreshDefault =  0.6;
pupilRangeDefault = [30 90];
pupilEllipseThreshDefault = 0.95;
p.addParameter('frameRate', frameRateDefault, @isnumeric);
p.addParameter('gammaCorrection', gammaCorrectionDefault, @isnumeric);
p.addParameter('pupilCircleThresh', pupilCircleThreshDefault, @isnumeric);
p.addParameter('pupilRange', pupilRangeDefault, @isnumeric);
p.addParameter('pupilEllipseThresh', pupilEllipseThreshDefault, @isnumeric);

%parse
p.parse(grayI, perimeterVideo, varargin{:})

% define optional variables values
frameRate = p.Results.frameRate;
gammaCorrection = p.Results.gammaCorrection;
pupilCircleThresh =  p.Results.pupilCircleThresh;
pupilRange = p.Results.pupilRange;
pupilEllipseThresh = p.Results.pupilEllipseThresh;


%% initiate output video object

outObj              = VideoWriter(perimeterVideoPath);
outObj.FrameRate    = frameRate;
open(outObj);

%% extract pupil perimeter

% get number of frames from grayI
numFrames = size(grayI,3);

% TAKING OUT PROGRESS BAR FOR NOW
% % initiate progress bar 
% progBar = ProgressBar(numFrames,'tracking pupil...');

% open a figure
ih = figure;

% loop through gray frames
for i = 1:numFrames
    
    % Get the frame
    I = squeeze(grayI(:,:,i));
    
     % adjust gamma for this frame
    I = imadjust(I,[],[],gammaCorrection);
    
    % Show the frame
    imshow(I);
    
    % track with circles 
    % feed default values for glint (glint data won't be stored now)
    glintCircleThresh = 0.999;
    glintRange = [10 30];
    
    [pCenters, pRadii,pMetric,~,~,~, pupilRange, ~] = circleFit(I,pupilCircleThresh,glintCircleThresh,pupilRange,glintRange,'pupilOnly', true);
    
    if isempty(pCenters) %no pupil circle patch was found
        % just save frame
        frame   = getframe(ih);
        writeVideo(outObj,frame);
%         % increment progress bar
%         if ~mod(i,10);progBar(i);end;
        continue
    else
        % get pupil perimeter
        [binP] = getPupilPerimeter(I,pCenters,pRadii, pupilEllipseThresh);
        
        % show the thresholded frame
        imshow(binP)
        
        % write frame to thresholded video
        frame   = getframe(ih);
        writeVideo(outObj,frame);
        
%         % increment progress bar
%         if ~mod(i,10);progBar(i);end;
    end
end % loop through gray frames

%% close video
close(ih);
close(outObj);

%% export perimeter params
perimeterParams = p.Results;
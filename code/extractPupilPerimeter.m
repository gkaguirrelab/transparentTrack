function [perimeter] = extractPupilPerimeter(grayVideoName, perimeterVideoName,varargin)
% function [perimeterParams] = extractPupilPerimeter(grayVideoName, perimeterVideoName,varargin)
% 
% This function thresholds the video to extract the pupil perimeter and
% saves out a BW video showing the pupil perimeter only.
% 
% Output
%   perimeter - matfile with informations about the parameters and the
%       environment in which the code was run.
% 
% Input (required)
%       grayVideoName - full path to  the gray video to track
%       perimeterVideoName - full path to the output avi in which to save 
%           the output.
%       
% Options (analysis)
%       gammaCorrection - gamma correction to be applied in current frame 
%           (default 1, typical range [0.5 1.8])
%       pupilCircleThresh - threshold value to locate the glint for circle
%           fitting (default 0.06, typical range [0.04 0.09])
%       pupilRange - initial radius range for circle fitting of the glint
%       (default [30 90]). This value gets dynamically updated.
%       pupilEllipseThresh - threshold value to locate the glint for
%           ellipse fitting (default 0.9, typical range [0.8 0.98])
%       nFrames - number of frames to process. If not specified or
%           Inf will process the full video.
% 
% Options (environment)
%   tbSnapshot - the passed tbSnapshot output that is to be saved along
%      with the data
%   timestamp / username / hostname - these are automatically derived and
%      saved within the p.Results structure. 
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
p.addParameter('gammaCorrection', 1, @isnumeric);
p.addParameter('pupilCircleThresh', 0.06, @isnumeric);
p.addParameter('pupilRange', [30 90], @isnumeric);
p.addParameter('pupilEllipseThresh', 0.95, @isnumeric);
p.addParameter('glintCircleThresh', 0.999, @isnumeric);
p.addParameter('glintRange', [10 30], @isnumeric);
p.addParameter('nFrames',Inf,@isnumeric);

% Optional display and I/O params
p.addParameter('verbosity','none',@ischar);
p.addParameter('showTracking', false, @islogical)

% Environment parameters
p.addParameter('tbSnapshot',[],@(x)(isempty(x) | isstruct(x)));
p.addParameter('timestamp',char(datetime('now')),@ischar);
p.addParameter('username',char(java.lang.System.getProperty('user.name')),@ischar);
p.addParameter('hostname',char(java.net.InetAddress.getLocalHost.getHostName),@ischar);

% parse
p.parse(grayVideoName, perimeterVideoName, varargin{:})

%% read input video file
% load pupilPerimeter
inObj = VideoReader(grayVideoName);

% get number of frames
if p.Results.nFrames == Inf
    nFrames = floor(inObj.Duration*inObj.FrameRate);
else
    nFrames = p.Results.nFrames;
end
%% initiate output video object

outObj = VideoWriter(perimeterVideoName);
outObj.FrameRate = inObj.FrameRate;
open(outObj);

%% extract pupil perimeter

% alert the user
if strcmp(p.Results.verbosity,'full')
    tic
    fprintf(['Extracting pupil perimeter. Started ' char(datetime('now')) '\n']);
    fprintf('| 0                      50                   100%% |\n');
    fprintf('.');
end

if p.Results.showTracking
% open a figure
ih = figure;
end

% initialize pupil range (it will change dynamically in the loop)
pupilRange = p.Results.pupilRange;

% loop through gray frames
for ii = 1:nFrames
    % increment the progress bar
    if strcmp(p.Results.verbosity,'full') && mod(ii,round(nFrames/50))==0
        fprintf('.');
    end
    
    % Get the frame
    I = readFrame(inObj);
    
     % adjust gamma for this frame
    I = imadjust(I,[],[],p.Results.gammaCorrection);
    
    I = rgb2gray (I);
    
    if p.Results.showTracking
    % show the frame
    imshow(I, 'Border', 'tight');
    end
    
    % track with circles 
    [pCenters, pRadii,~,~,~,~, pupilRange, ~] = circleFit(I,p.Results.pupilCircleThresh,p.Results.glintCircleThresh,pupilRange,p.Results.glintCircleThresh);
    
    if isempty(pCenters) %no pupil circle patch was found
        % make the frame black and save it
        I = zeros(size(I));
        
        if p.Results.showTracking
            imshow(I, 'Border', 'tight');
            frame   = getframe(ih);
        else
            frame   = im2uint8(I);
        end
        
        writeVideo(outObj,frame);
        continue
    else
        % get pupil perimeter
        [binP] = getPupilPerimeter(I,pCenters,pRadii, p.Results.pupilEllipseThresh);
        
        % convert BW frame
        I = im2uint8(binP);
        
        if p.Results.showTracking
            imshow(I, 'Border', 'tight');
            frame   = getframe(ih);
        else
            frame   = im2uint8(I);
        end

        writeVideo(outObj,frame);
        
    end
end % loop through gray frames

%% close video
if p.Results.showTracking
close(ih);
end
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

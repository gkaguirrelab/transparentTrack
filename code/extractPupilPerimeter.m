function perimeterParams = extractPupilPerimeter(grayI, perimeterVideoName,varargin)

% This function thresholds the video to extract the pupil perimeter and
% saves out a BW video showing the pupil perimeter only.
% 
% 
% Input params
% ============
%       grayI : 3D array of gray frames to track
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
%  define options beforehand:
% 
%  gammaCorrection = 1.2;
%  pupilEllipseThresh = 0.93;
%  perimeterParams = extractPupilPerimeter(grayI, perimeterVideoName, 'gammaCorrection', gammaCorrection, 'pupilEllipseThresh', pupilEllipseThresh)
% 
%  OR
% 
%  define options directly:
% 
%  perimeterParams = extractPupilPerimeter(grayI, perimeterVideoName, 'gammaCorrection', 1.2, 'pupilEllipseThresh', 0.93)


%% parse input and define variables

p = inputParser;
% required input
p.addRequired('grayI');
p.addRequired('perimeterVideo',@isstr);

% optional inputs
p.addParameter('frameRate', 60, @isnumeric);
p.addParameter('gammaCorrection', 1, @isnumeric);
p.addParameter('pupilCircleThresh', 0.06, @isnumeric);
p.addParameter('pupilRange', [30 90], @isnumeric);
p.addParameter('pupilEllipseThresh', 0.95, @isnumeric);
p.addParameter('verbosity','none', @ischar);

% parse
p.parse(grayI, perimeterVideoName, varargin{:})

% define optional variables values
frameRate = p.Results.frameRate;
gammaCorrection = p.Results.gammaCorrection;
pupilCircleThresh =  p.Results.pupilCircleThresh;
pupilRange = p.Results.pupilRange;
pupilEllipseThresh = p.Results.pupilEllipseThresh;


%% initiate output video object
outObj = VideoWriter(perimeterVideoName);
outObj.FrameRate = frameRate;
open(outObj);

%% extract pupil perimeter

% get number of frames from grayI
nFrames = size(grayI,3);

% alert the user
if strcmp(p.Results.verbosity,'full')
    tic
    fprintf(['Extracting pupil perimeter. Started ' char(datetime('now')) '\n']);
    fprintf('| 0                      50                   100%% |\n');
    fprintf('.');
end

% open a figure
ih = figure;

% loop through gray frames
for ii = 1:nFrames
    % increment the progress bar
    if strcmp(p.Results.verbosity,'full') && mod(ii,round(nFrames/50))==0
        fprintf('.');
    end
    
    % Get the frame
    I = squeeze(grayI(:,:,ii));
    
     % adjust gamma for this frame
    I = imadjust(I,[],[],gammaCorrection);
    
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
        
    end
end % loop through gray frames

%% close video
close(ih);
clear outObj

%% export perimeter params
perimeterParams = p.Results;

% report completion of analysis
if strcmp(p.Results.verbosity,'full')
    fprintf('\n');
    toc
    fprintf('\n');
end

end % function
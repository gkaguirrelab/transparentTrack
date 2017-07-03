function [perimeter] = extractPupilPerimeter(grayVideoName, perimeterFileName, varargin)
% function [perimeterParams] = extractPupilPerimeter(grayVideoName, perimeterVideoName,varargin)
% 
% This function thresholds the video to extract the pupil perimeter and
% saves out a BW video showing the pupil perimeter only.
%
% An initial search for the pupil border is performed with the circleFit
% function. If a candidate circle is found, the region is dilated. We then
% binarize the resulting "patch" image with a user determined threshold,
% and extract the perimeter of the bigger region surviving the thresholding
% process (believed to be the pupil).
% 
% Output
%   perimeter - structure with a 'data' field containing a 3D matrix
%       (video height x video width x nFrames) containing the pupil
%       perimeter, and a meta field with analysis and environment params.
% 
% Input (required)
%	grayVideoName - full path to  the gray video to track
%	perimeterFileName - full path to the .mat file in which to save the
%       output.
%       
% Options (analysis)
% 	gammaCorrection - gamma correction to be applied to the video frames 
%       (default 1, typical range [0.5 1.8])
%   pupilCircleThresh - threshold value to locate the glint for circle
%       fitting (default 0.06, typical range [0.04 0.09])
%	pupilRange - initial radius range for circle fitting of the glint
%       (default [30 90]). This value gets dynamically updated.
%	pupilEllipseThresh - threshold value to locate the glint for
%       ellipse fitting (default 0.9, typical range [0.8 0.98])
%   glintCircleThresh - DEFINE HERE
%   glintRange - DEFINE HERE
%   maskBox - DEFINE HERE
%   smallObjThresh - DEFINE HERE
%   nFrames - number of frames to process. If not specified or
%       Inf will process the full video.
% 
% Options (display)
%   verbosity - controls console status updates
%   showTracking - controls display of the current frame and perimeter
%
% Options (environment)
%   tbSnapshot - the passed tbSnapshot output that is to be saved along
%      with the data
%   timestamp / username / hostname - these are automatically derived and
%      saved within the p.Results structure. 
% 


%% Parse input and define variables
p = inputParser;

% required input
p.addRequired('grayVideoName',@isstr);
p.addRequired('perimeterFileName',@isstr);

% Optional analysis params
p.addParameter('gammaCorrection', 1, @isnumeric);
p.addParameter('pupilCircleThresh', 0.06, @isnumeric);
p.addParameter('pupilRange', [30 90], @isnumeric);
p.addParameter('pupilEllipseThresh', 0.95, @isnumeric);
p.addParameter('glintCircleThresh', 0.999, @isnumeric);
p.addParameter('glintRange', [10 30], @isnumeric);
p.addParameter('maskBox', [4 30], @isnumeric);
p.addParameter('smallObjThresh', 500, @isnumeric);
p.addParameter('nFrames',Inf,@isnumeric);

% Optional display params
p.addParameter('verbosity','none',@ischar);
p.addParameter('showTracking', false, @islogical)

% Environment parameters
p.addParameter('tbSnapshot',[],@(x)(isempty(x) | isstruct(x)));
p.addParameter('timestamp',char(datetime('now')),@ischar);
p.addParameter('username',char(java.lang.System.getProperty('user.name')),@ischar);
p.addParameter('hostname',char(java.net.InetAddress.getLocalHost.getHostName),@ischar);

% parse
p.parse(grayVideoName, perimeterFileName, varargin{:})


%% Display setup and prepare video objects
% alert the user
if strcmp(p.Results.verbosity,'full')
    tic
    fprintf(['Extracting pupil perimeter. Started ' char(datetime('now')) '\n']);
    fprintf('| 0                      50                   100%% |\n');
    fprintf('.');
end

% open a figure
if p.Results.showTracking
    figureHandle = figure;
end

% open inVideoObject
inVideoObj = VideoReader(grayVideoName);
% get number of frames
if p.Results.nFrames == Inf
    nFrames = floor(inVideoObj.Duration*inVideoObj.FrameRate);
else
    nFrames = p.Results.nFrames;
end

% get video dimensions
videoSizeX = inVideoObj.Width;
videoSizeY = inVideoObj.Height;

% initialize variable to hold the perimeter data
perimeter.data = zeros(videoSizeY,videoSizeX,nFrames,'uint8');


%% Extract pupil perimeter
% initialize pupil range (it will change dynamically in the loop)
pupilRange = p.Results.pupilRange;

% structuring element for pupil mask size
sep = strel('rectangle',p.Results.maskBox);

% loop through gray frames
for ii = 1:nFrames

    % increment the progress bar
    if strcmp(p.Results.verbosity,'full') && mod(ii,round(nFrames/50))==0
        fprintf('.');
    end
    
    % Read the frame, adjust gamma, make gray
    thisFrame = readFrame(inVideoObj);
    thisFrame = imadjust(thisFrame,[],[],p.Results.gammaCorrection);
    thisFrame = rgb2gray (thisFrame);

    % show the initial frame
    if p.Results.showTracking
        imshow(thisFrame, 'Border', 'tight');
    end
    
    % perform an initial search for the pupil with circleFit 
    [pCenters, pRadii,~,~,~,~, pupilRange, ~] = ...
        circleFit(thisFrame,...
        p.Results.pupilCircleThresh,...
        p.Results.glintCircleThresh,...
        pupilRange,...
        p.Results.glintCircleThresh);
    
    % If a pupil circle patch was found, get the perimeter, else write out
    % a zero-filled frame
    if ~isempty(pCenters)
        
        % generate mask
        pupilMask = zeros(size(thisFrame));
        pupilMask = insertShape(pupilMask,'FilledCircle',[pCenters(1,1) pCenters(1,2) pRadii(1)],'Color','white');
        pupilMask = imdilate(pupilMask,sep);
        pupilMask = im2bw(pupilMask);
        
        % apply mask to grey image complement image
        complementThisFrame = imcomplement(thisFrame);
        maskedPupil = immultiply(complementThisFrame,pupilMask);
        
        % convert back to gray
        pI = uint8(maskedPupil);
        
        % Binarize pupil
        binP = ones(size(pI));
        binP(pI<quantile(double(complementThisFrame(:)),p.Results.pupilEllipseThresh)) = 0;
        
        % remove small objects
        binP = bwareaopen(binP, p.Results.smallObjThresh);
        
        % fill the holes
        binP = imfill(binP,'holes');
        
        % get perimeter of object
        binP = bwperim(binP);
        
        % save the perimeter
        thisFrame = im2uint8(binP);
        perimeter.data(:,:,ii) = thisFrame;
    else
        thisFrame = im2uint8(zeros(size(thisFrame)));
        perimeter.data(:,:,ii) = thisFrame;
    end
    
    % show the perimeter frame
    if p.Results.showTracking
        imshow(thisFrame, 'Border', 'tight');
    end
        
end % loop through gray frames

%% Clean up and close
if p.Results.showTracking
    close(figureHandle);
end

clear inVideoObj
clear outVideoObj

% save mat file with the video and analysis details
perimeter.meta = p.Results;
save(perimeterFileName,'perimeter');

% report completion of analysis
if strcmp(p.Results.verbosity,'full')
    fprintf('\n');
    toc
    fprintf('\n');
end

end % function

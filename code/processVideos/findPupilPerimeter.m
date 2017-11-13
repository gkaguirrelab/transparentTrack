function perimeter = findPupilPerimeter(grayVideoName, perimeterFileName, varargin)
% function findPupilPerimeter(grayVideoName, perimeterVideoName,varargin)
%
% This function thresholds the video to find the pupil perimeter.
%
% An initial search for the pupil border is performed with the local
% function 'findPupilCircle'. If a candidate circle is found, the region is
% dilated. We then binarize the resulting "patch" image with a user
% set threshold, and extract the perimeter of the bigger region surviving 
% the thresholding process (believed to be the pupil).
%
% Output
%   perimeter - structure with a 'data' field containing a 3D matrix
%       (video height x video width x nFrames) containing the pupil
%       perimeter, and a meta field with analysis and environment params.
%
% Input (required)
%	grayVideoName - full path to the gray video to track
%	perimeterFileName - full path to the .mat file in which to save the
%       output.
%
% Options (analysis)
% 	pupilGammaCorrection - gamma correction to be applied to the video frames
%       (default 1, typical range [0.5 1.8])
%   pupilFrameMask : this option with add a mask on the original gray video,
%       framing it by [nRows nColumns] on the borders symmetrically or by
%       [nRowsTop nColumnsRight nRowsBottom nColumnsLeft]. This is
%       particularly useful for size calibration videos in which appear
%       partial black dots that may throw off the circle finding mechanism.
%       fitting (default 0.06, typical range [0.04 0.09])
%   maskBox - This is the proportion to dilate the pupil masked region in
%       the vertical and horizontal directions respectively. A value of
%       zero will result in no dilation in that direction. A value of unity
%       will result in a masked region that is twice the size of the pupil
%       radius.
%   frameMaskValue - the image value that is assigned to the region that is
%       masked by frameMask. This should be a gray that is neither pupil
%       nor glint.
%   smallObjThresh - maximum size of small objects to be removed to clean
%       up the pupil perimeter.
%
% Options (used in the local function findPupilCircle)
%   pupilCircleThresh - The threshold used to binarize the image. Ranges
%      from 0 to 1, with 0 corresponding to black.
%	pupilRange - initial radius range in pixels for circle fitting of the
%       pupil (default [20 189]). This value gets dynamically updated.
%   imfindcirclesSensitivity - a parameter (ranging from 0-1) that is used
%       in the call to the matlab function imfindcircles as part of the
%       key-value pair 'Sensitivity'. We set the default value to 0.99,
%       meaning that the function will be very sensitive to circles, even
%       those that are partially obscured.
%   rangeAdjust - for each frame, a circle is fit to the pupil. The size of
%       the found pupil is used to update the pupilRange values for the
%       next frame. The rangeAdjust parameter defines the proportion above
%       and below this target size that is used to search on the next frame
%
% Options (verbosity and display)
%   verbosity - controls console status updates
%   displayMode - when set to true, displays the results of the boundary
%   	extraction and does not save a video.
%
% Optional key/value pairs (flow control)
%  'nFrames' - analyze fewer than the total number of frames.
%  'startFrame' - which frame to start on
%
% Options (environment)
%   tbSnapshot - the passed tbSnapshot output that is to be saved along
%      with the data
%   timestamp / username / hostname - these are automatically derived and
%      saved within the p.Results structure.
%


%% Parse input and define variables
p = inputParser; p.KeepUnmatched = true;

% required input
p.addRequired('grayVideoName',@isstr);
p.addRequired('perimeterFileName',@isstr);

% Optional analysis params
p.addParameter('pupilGammaCorrection', 0.75, @isnumeric);
p.addParameter('pupilFrameMask', [], @isnumeric);
p.addParameter('maskBox', [0.20 0.75], @isnumeric);
p.addParameter('frameMaskValue', 220, @isnumeric);
p.addParameter('smallObjThresh', 400, @isnumeric);

% findPupilCircle routine params. Defined here for transparency
p.addParameter('pupilCircleThresh', 0.06, @isnumeric);
p.addParameter('pupilRange', [40 800], @isnumeric);
p.addParameter('imfindcirclesSensitivity', 0.99, @isnumeric);
p.addParameter('rangeAdjust', 0.05, @isnumeric);

% Optional display params
p.addParameter('verbosity','none',@ischar);
p.addParameter('displayMode',false,@islogical);

% Optional flow control params
p.addParameter('nFrames',Inf,@isnumeric);
p.addParameter('startFrame',1,@isnumeric);

% Environment parameters
p.addParameter('tbSnapshot',[],@(x)(isempty(x) | isstruct(x)));
p.addParameter('timestamp',char(datetime('now')),@ischar);
p.addParameter('username',char(java.lang.System.getProperty('user.name')),@ischar);
p.addParameter('hostname',char(java.net.InetAddress.getLocalHost.getHostName),@ischar);

% parse
p.parse(grayVideoName, perimeterFileName, varargin{:})


%% Read video file into memory
% load pupilPerimeter
videoInObj = VideoReader(grayVideoName);
% get number of frames
if p.Results.nFrames == Inf
    nFrames = floor(videoInObj.Duration*videoInObj.FrameRate);
else
    nFrames = p.Results.nFrames;
end
% get video dimensions
videoSizeX = videoInObj.Width;
videoSizeY = videoInObj.Height;
% initialize variable to hold the perimeter data
grayVideo = zeros(videoSizeY,videoSizeX,nFrames,'uint8');
% read the video into memory, adjusting gamma and local contrast
for ii = 1:nFrames
    thisFrame = readFrame(videoInObj);
    thisFrame = imadjust(thisFrame,[],[],p.Results.pupilGammaCorrection);
    grayVideo(:,:,ii) = rgb2gray (thisFrame);
end
% close the video object
clear videoInObj


%% Extract pupil perimeter

% Detect display mode
if p.Results.displayMode
    fprintf('** DISPLAY MODE **\n')
    fprintf('Results will not be saved. Press space at any time to quit routine.\n')
    
    % create a figure for display
    figureHandle=figure();
    
    % we will monitor the currentchar for a 'q'
    set(figureHandle,'currentchar','?')
end

% alert the user
if strcmp(p.Results.verbosity,'full')
    tic
    fprintf(['Extracting pupil perimeter. Started ' char(datetime('now')) '\n']);
    fprintf('| 0                      50                   100%% |\n');
    fprintf('.');
end

% initialize variable to hold the perimeter data
perimeter_data = zeros(videoSizeY,videoSizeX,nFrames,'uint8');

% Initialize the pupilRange with the parameter value. This value is updated
% as we progress through the frames
pupilRange= p.Results.pupilRange;

% loop through gray frames
for ii = p.Results.startFrame:nFrames
    
    if p.Results.displayMode && strcmp(get(figureHandle,'currentchar'),' ')
        close(figureHandle)
        return
    end
    
    % increment the progress bar
    if strcmp(p.Results.verbosity,'full') && mod(ii-p.Results.startFrame+1,round(nFrames/50))==0
        fprintf('.');
    end
    
    % get the frame
    thisFrame = squeeze(grayVideo(:,:,ii));
    
   % apply a frame mask if required
    if ~isempty (p.Results.pupilFrameMask)
        if length(p.Results.pupilFrameMask) == 2
            thisFrame((1:p.Results.pupilFrameMask(1)),:) = p.Results.frameMaskValue;
            thisFrame((end - p.Results.pupilFrameMask(1):end),:) = p.Results.frameMaskValue;
            thisFrame(:, (1:p.Results.pupilFrameMask(2))) = p.Results.frameMaskValue;
            thisFrame(:, (end - p.Results.pupilFrameMask(2):end)) = p.Results.frameMaskValue;
        elseif length(p.Results.pupilFrameMask) == 4
            thisFrame((1:p.Results.pupilFrameMask(1)),:) = p.Results.frameMaskValue; %top
            thisFrame(:, (end - p.Results.pupilFrameMask(2):end)) = p.Results.frameMaskValue; %left
            thisFrame((end - p.Results.pupilFrameMask(3):end),:) = p.Results.frameMaskValue; %bottom
            thisFrame(:, (1:p.Results.pupilFrameMask(4))) = p.Results.frameMaskValue; %right
        else
            error ('invalid frameMask parameter. Frame mask must be defined as [nRows nColumns] or as [nRowsTop nColumnsRight nRowsBottom nColumnsLeft]')
        end
    end
    
    % store the current pupilRange
    initialPupilRange = pupilRange;    
    
    % perform an initial search for the pupil with findGlintAndPupilCircles. Also extract
    % glint location and size information for later use.
    [pCenters, pRadii,~,pupilRange] = ...
        findPupilCircle(thisFrame,...
        p.Results.pupilCircleThresh,...
        pupilRange,...
        p.Results.imfindcirclesSensitivity,p.Results.rangeAdjust);
    
    % If a pupile circle patch was not found, try again after expanding the
    % pupil search range by 50%, then 100%. We limit the possible range for
    % the pupil search to the passed default bounds
    if isempty(pCenters)
        % Check if a 50% expansion is within bounds
        candidateRange = [ceil(initialPupilRange(1)/1.5) round(initialPupilRange(2)*1.5)];
        if candidateRange(1) < p.Results.pupilRange(1) || candidateRange(2) > p.Results.pupilRange(2)
            candidateRange = p.Results.pupilRange;
        end
        [pCenters, pRadii,~,pupilRange] = ...
            findPupilCircle(thisFrame,...
            p.Results.pupilCircleThresh,...
            candidateRange,...
            p.Results.imfindcirclesSensitivity,p.Results.rangeAdjust);
        if isempty(pCenters) % still no circle? Try 100% increase
            candidateRange = [ceil(initialPupilRange(1)/2) round(initialPupilRange(2)*2)];
            if candidateRange(1) < p.Results.pupilRange(1) || candidateRange(2) > p.Results.pupilRange(2)
                candidateRange = p.Results.pupilRange;
            end
            [pCenters, pRadii,~,pupilRange] = ...
                findPupilCircle(thisFrame,...
                p.Results.pupilCircleThresh,...
                candidateRange,...
                p.Results.imfindcirclesSensitivity,p.Results.rangeAdjust);
            if isempty(pCenters) % STILL no circle? Give up and restore initialPupilRange
                pupilRange = initialPupilRange;
            end
        end
    end
    
    % If a pupil circle patch was ultimately found, get the perimeter, else
    % write out a zero-filled frame
    if ~isempty(pCenters)
        
        % structuring element for pupil mask size. This is a rectangular
        % dilation box that is adapted to the size of the radius of the
        % initially found pupil circle. The proportional size of the
        % dilation box is set in the key value 'maskBox'.
        sep = strel('rectangle',round(pRadii(1).*p.Results.maskBox));
        
        % generate mask
        pupilMask = zeros(size(thisFrame));
        pupilMask = insertShape(pupilMask,'FilledCircle',[pCenters(1,1) pCenters(1,2) pRadii(1)],'Color','white');
        pupilMask = imdilate(pupilMask,sep);
        pupilMask = im2bw(pupilMask);
        
        % apply mask to the complement of the image
        complementThisFrame = imcomplement(thisFrame);
        maskedPupil = immultiply(complementThisFrame,pupilMask);
        
        % partition the image into three regions, corresponding (from
        % lightest to darkest) to the glint, iris, and pupil. We first NaN
        % out all points in the image that are not within the masked
        % region, so that they do not influence the three region partition.
        maskedPupilNaN=double(maskedPupil);
        maskedPupilNaN(pupilMask==0)=NaN;
        
        % save the warning state, then turn off the warning that can occur
        % if multiThresh is unable to find three regions
        warningState = warning;
        warning('off','images:multithresh:degenerateInput');
        
        % perform the multi thresh
        otsuThresh = multithresh(maskedPupilNaN,2);
        
        % restore warning state
        warning(warningState);
                
        % Set the glint and iris to zero, and the pupil to unity
        binP = imquantize(maskedPupil, otsuThresh, [0 0 1]);
        
        % remove small objects
        binP = bwareaopen(binP, p.Results.smallObjThresh);
        
        % fill the holes
        binP = imfill(binP,'holes');
        
        % get perimeter of object
        binP = bwperim(binP);
        
        % save the perimeter
        perimFrame = im2uint8(binP);
        perimeter_data(:,:,ii) = perimFrame;
    else
        perimFrame = im2uint8(zeros(size(thisFrame)));
        perimeter_data(:,:,ii) = perimFrame;
        pupilRange = initialPupilRange;
    end
    
    if p.Results.displayMode
        displayFrame=thisFrame;
        [Yp, Xp] = ind2sub(size(perimFrame),find(perimFrame));
        if ~isempty(Xp)
            displayFrame(sub2ind(size(perimFrame),Yp,Xp))=255;
        end
        imshow(displayFrame, 'Border', 'tight');
    end
    
end % loop through gray frames

%% Clean up and close

% save mat file with the video and analysis details
perimeter.data = perimeter_data;
perimeter.meta = p.Results;
if ~p.Results.displayMode
    save(perimeterFileName,'perimeter','-v7.3');
else
    close(figureHandle);
end

% report completion of analysis
if strcmp(p.Results.verbosity,'full')
    fprintf('\n');
    toc
    fprintf('\n');
end


end %  main function


function [pCenters, pRadii,pMetric, pupilRange] = findPupilCircle(I,pupilCircleThresh,pupilRange,imfindcirclesSensitivity,rangeAdjust)
% findGlintAndPupilCircles(I,pupilCircleThresh,glintCircleThresh,pupilRange,glintRange,pupilOnly,glintOut,dilateGlint,imfindcirclesSensitivity,rangeAdjust)
%
% this function is used for pupil circle fitting.

%% parse input and define variables
p = inputParser;
% required input
p.addRequired('I');
p.addRequired('pupilCircleThresh', @isnumeric)
p.addRequired('pupilRange',@isnumeric);
p.addRequired('imfindcirclesSensitivity', @isnumeric);
p.addRequired('rangeAdjust', @isnumeric);

% parse
p.parse(I,pupilCircleThresh,pupilRange,imfindcirclesSensitivity,rangeAdjust);

%% circle fit

% create blurring filter
filtSize = round([0.01*min(size(I)) 0.01*min(size(I)) 0.01*min(size(I))]);


% Filter for pupil
padP = padarray(I,[size(I,1)/2 size(I,2)/2], 128);
h = fspecial('gaussian',[filtSize(1) filtSize(2)],filtSize(3));
pI = imfilter(padP,h);
pI = pI(size(I,1)/2+1:size(I,1)/2+size(I,1),size(I,2)/2+1:size(I,2)/2+size(I,2));
% Binarize pupil
binP = ones(size(pI));
binP(pI<quantile(double(pI(:)),pupilCircleThresh)) = 0;


% store the warning state
origWarnState = warning;

% Silence the imfindcircles warning regarding circle size
warning('off','images:imfindcircles:warnForLargeRadiusRange');
warning('off','images:imfindcircles:warnForSmallRadius');

% Find the pupil
[pCenters, pRadii,pMetric] = imfindcircles(binP,pupilRange,'ObjectPolarity','dark',...
    'Sensitivity',imfindcirclesSensitivity);

% Restore the warning state
warning(origWarnState);


% adjust the pupil range (for quicker processing)
if ~isempty(pCenters)
    pupilRange(1)   = min(floor(pRadii(1)*(1-rangeAdjust)),pupilRange(2));  
    pupilRange(2)   = max(ceil(pRadii(1)*(1 + rangeAdjust)),pupilRange(1));
else
    pupilRange(1)   = max(ceil(pupilRange(1)*(1 - rangeAdjust)),pupilRange(1));
    pupilRange(2)   = min(ceil(pupilRange(2)*(1 + rangeAdjust)),pupilRange(2));
end

end % function

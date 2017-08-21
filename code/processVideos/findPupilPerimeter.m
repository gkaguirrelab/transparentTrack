function findPupilPerimeter(grayVideoName, perimeterFileName, varargin)
% function findPupilPerimeter(grayVideoName, perimeterVideoName,varargin)
%
% This function thresholds the video to find the pupil perimeter.
%
% An initial search for the pupil border is performed with the
%
%       findPupilCircle
%
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
% 	pupilGammaCorrection - gamma correction to be applied to the video frames
%       (default 1, typical range [0.5 1.8])
%   pupilCircleThresh - threshold value to locate the pupil for circle
%       fitting (default 0.06, typical range [0.04 0.09])
%	pupilRange - initial radius range for circle fitting of the pupil
%       (default [30 90]). This value gets dynamically updated.
%   maskBox - This is the proportion to dilate the pupil masked region in
%       the vertical and horizontal directions respectively. A value of
%       zero will result in no dilation in that direction. A value of unity
%       will result in a masked region that is twice the size of the pupil
%       radius.
%   frameMask - this option with add a mask on the original gray video, 
%       framing it by [nRows nColumns] on the borders. This is particularly 
%       useful for size calibration videos in which appear partial black 
%       dots that may throw off the circle finding mechanism.
%   frameMaskValue - the image value that is assigned to the region that is
%       masked by frameMask. This should be a gray that is neither pupil
%       nor glint.
%   smallObjThresh - maximum size of small objects to be removed to clean
%       up the pupil perimeter.
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
p.addParameter('maskBox', [0.20 0.75], @isnumeric);
p.addParameter('frameMask', [], @isnumeric);
p.addParameter('frameMaskValue', 220, @isnumeric);
p.addParameter('smallObjThresh', 200, @isnumeric);

% findPupilCircle routine params. Defined here for transparency
p.addParameter('pupilCircleThresh', 0.06, @isnumeric);
p.addParameter('pupilRange', [20 180], @isnumeric);
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
    
    % apply the frameMaskValue to the frame region (if so requested)
    if ~isempty (p.Results.frameMask)
        thisFrame((1:p.Results.frameMask(1)),:) = p.Results.frameMaskValue;
        thisFrame((end - p.Results.frameMask(1):end),:) = p.Results.frameMaskValue;
        thisFrame(:, (1:p.Results.frameMask(2))) = p.Results.frameMaskValue;
        thisFrame(:, (end - p.Results.frameMask(2):end)) = p.Results.frameMaskValue;
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
    save(perimeterFileName,'perimeter');
else
    close(figureHandle);
end

% report completion of analysis
if strcmp(p.Results.verbosity,'full')
    fprintf('\n');
    toc
    fprintf('\n');
end


end % function

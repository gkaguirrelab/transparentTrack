function extractPupilPerimeter(grayVideoName, perimeterFileName, varargin)
% function extractPupilPerimeter(grayVideoName, perimeterVideoName,varargin)
% 
% This function thresholds the video to extract the pupil perimeter.
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
%   glintCircleThresh - DEFINE HERE
%   glintRange - DEFINE HERE
%   maskBox - This is the proportion to dilate the pupil masked region in
%       the vertical and horizontal directions respectively. A value of
%       zero will result in no dilation in that direction. A value of unity
%       will result in a masked region that is twice the size of the pupil
%       radius.
%   smallObjThresh - DEFINE HERE
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
p.addParameter('gammaCorrection', 1, @isnumeric);
p.addParameter('pupilCircleThresh', 0.06, @isnumeric);
p.addParameter('pupilRange', [30 90], @isnumeric);
p.addParameter('glintCircleThresh', 0.999, @isnumeric);
p.addParameter('glintRange', [10 30], @isnumeric);
p.addParameter('maskBox', [0.20 0.75], @isnumeric);
p.addParameter('smallObjThresh', 500, @isnumeric);
p.addParameter('adaptHisEq', true, @islogical);
p.addParameter('localContrastEdgeThresh', 0.5, @isnumeric);
p.addParameter('localContrastAmount', 0.5, @isnumeric);

% circleFit routine params. Defined here for transparency
p.addParameter('pupilOnly', false, @islogical);
p.addParameter('glintOut', 0.1, @isnumeric);
p.addParameter('dilateGlint', 6, @isnumeric);
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
    thisFrame = imadjust(thisFrame,[],[],p.Results.gammaCorrection);
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
            
    % perform an initial search for the pupil with circleFit 
    [pCenters, pRadii,~,~,~,~, pupilRange, ~] = ...
        circleFit(thisFrame,...
        p.Results.pupilCircleThresh,...
        p.Results.glintCircleThresh,...
        pupilRange,...
        p.Results.glintRange,...
        p.Results.pupilOnly,p.Results.glintOut,p.Results.dilateGlint,p.Results.imfindcirclesSensitivity,p.Results.rangeAdjust);
    
    % If a pupil circle patch was found, get the perimeter, else write out
    % a zero-filled frame
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
        
        % apply mask to grey image complement image
        complementThisFrame = imcomplement(thisFrame);
        maskedPupil = immultiply(complementThisFrame,pupilMask);

        % partition the image into three regions, corresponding (from
        % lightest to darkest) to the glint, iris, and pupil. We first NaN
        % out all points in the image that are not within the masked
        % region, so that they do not influence the three region partition.        
        maskedPupilNaN=double(maskedPupil);
        maskedPupilNaN(pupilMask==0)=NaN;
        otsuThresh = multithresh(maskedPupilNaN,2);
        
        % Now set the glint and iris to white and the pupil to black. 
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

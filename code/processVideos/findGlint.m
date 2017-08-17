function [glintData] = findGlint(grayVideoName, glintFileName, varargin)
%  [glintData] = trackGlint(grayVideoName, glintFileName, varargin)
%
% This function tracks one or more glints in the eye video using a simple
% thresholding and region property identification approach.
%
% Every frame is firstly corrected with an elevated gamma value, so that
% glints and other large bright spot are enhanced. The image is then
% binarized with a relatively high threshold, so that only the brights spot
% and their immediate surroundings results as non-zero values. The
% centroids of each surviving region in the binary image are extracted using
% matlab's function "regionpros". The centroid location is weighted with
% the actual brightness value of each pixel in the gray gamma-corrected
% image.
% 
% After all centroids location are extracted, data is refined according to
% the expected number of glints and average centroid location throughout
% the video. Firstly, we calculate the median location of the glints from
% those frames that return as many centroids as the desired glints. For
% those frames in which more than the expected number of glints is found,
% we use the median value of the "good centroids" location to assess which
% of the regions identified are indeed the desired glints. In frames where
% less than the desired number of glints is located, the missing centroids
% locations will be set as NaNs.
% 
% If the refineGlintLocation flag is set to true, the data is further
% cleaned up to exclude glint values that deviate from the median value
% more than the custom threshold maxGlintDeviation.
%
% Output
%	glintData - structure with fields that contain the X and Y location of
%       the center of the glint (in units of pixels), and a meta field with
%       additional analysis params and intermediate results.
%
% Input (required)
%	grayVideoName : full path to the video in which to track the glint.
%   	A grayscale video is expected.
%   glintFileName : full path to the matFile in which to save the glint
%       results.
%
% Options (analysis)
%   numberOfGlints : desired number of glint to find
%	glintGammaCorrection : gamma correction to be applied in current
%       frame. An extremely high value will make almost all the frame black
%       and only big bright spots will be white. This reduces the
%       possibility of confusing the glint with some other smaller bright
%       spot (default 4, decrease if no glint is found)
%   glintThreshold : threshold value to binarize the glint gray image.
%       Should be set to preserve both the glints and any "halo" around
%       them.(default value 0.8)
%   frameMask : 
%   centroidsAllocation : max number of centroids to be saved in memory

% Optional key/value pairs (flow control)
%  'nFrames' - analyze fewer than the total number of frames.
%  'useParallel' - If set to true, use the Matlab parallel pool for the
%    initial ellipse fitting.
%  'nWorkers' - Specify the number of workers in the parallel pool. If
%    undefined the default number will be used.
%  'tbtbProjectName' - The workers in the parallel pool are configured by
%    issuing a tbUseProject command for the project specified here.
%
% Options (environment)
%   tbSnapshot - the passed tbSnapshot output that is to be saved along
%      with the data
%   timestamp / username / hostname - these are automatically derived and
%      saved within the p.Results structure.
%

%% parse input and define variables

p = inputParser; p.KeepUnmatched = true;

% Required input
p.addRequired('grayVideoName',@isstr);
p.addRequired('glintFileName',@isstr);

% optional analysis parameters
p.addParameter('numberOfGlints', 1, @isnumeric);
p.addParameter('glintGammaCorrection', 4, @isnumeric);
p.addParameter('glintThreshold', 0.8, @isnumeric);
p.addParameter('frameMask',[] , @isnumeric);
p.addParameter('centroidsAllocation', 5, @isnumeric);

% Optional display params
p.addParameter('verbosity','none',@ischar);
p.addParameter('displayMode',false,@islogical);

% Optional flow control params
p.addParameter('nFrames',Inf,@isnumeric);
p.addParameter('useParallel',false,@islogical);
p.addParameter('nWorkers',[],@(x)(isempty(x) | isnumeric(x)));
p.addParameter('tbtbRepoName','LiveTrackAnalysisToolbox',@ischar);

% Environment parameters
p.addParameter('tbSnapshot',[],@(x)(isempty(x) | isstruct(x)));
p.addParameter('timestamp',char(datetime('now')),@ischar);
p.addParameter('username',char(java.lang.System.getProperty('user.name')),@ischar);
p.addParameter('hostname',char(java.net.InetAddress.getLocalHost.getHostName),@ischar);

% parse
p.parse(grayVideoName, glintFileName, varargin{:})


%% read video file into memory
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

% read the video into memory, adjust gamma
for ii = 1:nFrames
    thisFrame = readFrame(videoInObj);
    thisFrame = imadjust(thisFrame,[],[],p.Results.glintGammaCorrection);
    grayVideo(:,:,ii) = rgb2gray (thisFrame);
end
% close the video object
clear videoInObj

%% Set up the parallel pool
if p.Results.useParallel
    if strcmp(p.Results.verbosity,'full')
        tic
        fprintf(['Opening parallel pool. Started ' char(datetime('now')) '\n']);
    end
    if isempty(p.Results.nWorkers)
        parpool;
    else
        parpool(p.Results.nWorkers);
    end
    poolObj = gcp;
    if isempty(poolObj)
        nWorkers=0;
    else
        nWorkers = poolObj.NumWorkers;
        % Use TbTb to configure the workers.
        if ~isempty(p.Results.tbtbRepoName)
            spmd
                tbUse(p.Results.tbtbRepoName,'reset','full','verbose',false,'online',false);
            end
            if strcmp(p.Results.verbosity,'full')
                fprintf('CAUTION: Any TbTb messages from the workers will not be shown.\n');
            end
        end
    end
    if strcmp(p.Results.verbosity,'full')
        toc
        fprintf('\n');
    end
else
    nWorkers=0;
end

%% Initialize variables
% Initialize glint variables
glintData_X = nan(nFrames,p.Results.numberOfGlints);
glintData_Y = nan(nFrames,p.Results.numberOfGlints);

% initialize centroid variable according to the centroidsAllocation value
centroidsByFrame_X = nan(nFrames,p.Results.centroidsAllocation);
centroidsByFrame_Y = nan(nFrames,p.Results.centroidsAllocation);


%% Find all centroids

% Detect display mode
if p.Results.displayMode && ~p.Results.useParallel
    fprintf('** DISPLAY MODE **\n')
    fprintf('Results will not be saved.\n')
    
    % create a figure for display
    figureHandle=figure();
    
    % we will monitor the currentchar for a 'q'
    set(figureHandle,'currentchar','?')
end

% alert the user
if strcmp(p.Results.verbosity,'full')
    tic
    fprintf(['Tracking the glint. Started ' char(datetime('now')) '\n']);
    fprintf('| 0                      50                   100%% |\n');
    fprintf('.\n');
end


% loop through frames
% parfor (ii = 1:nFrames, nWorkers)
for   ii = 1:nFrames
    % increment the progress bar
    if strcmp(p.Results.verbosity,'full') && mod(ii,round(nFrames/50))==0
        fprintf('\b.\n');
    end
    
    % get the frame
    thisFrame = squeeze(grayVideo(:,:,ii));
    
    % apply a frame mask if required
    if ~isempty (p.Results.frameMask)
        thisFrame((1:p.Results.frameMask(1)),:) = 220;
        thisFrame((end - p.Results.frameMask(1):end),:) = 220;
        thisFrame(:, (1:p.Results.frameMask(2))) = 220;
        thisFrame(:, (end - p.Results.frameMask(2):end)) = 220;
    end
    
    % binarize glint image according to glintThreshold
    binG = im2bw(thisFrame, p.Results.glintThreshold);
    
    % get the weighted centroids for all the surviving bright spots in the
    % image. The "weight" is given by the actual gray value in the frame,
    % so that the brightest glint pixels are more relevant to the overall
    % glint position.
    stats = regionprops(binG, thisFrame, 'WeightedCentroid');
    
    % if centroids were found in this frame, save them out.
    if ~isempty(stats)
        centroids = stats.WeightedCentroid;
        for cc = 1: min(size(centroids,1),p.Results.centroidsAllocation)
            centroidsByFrame_X(ii,cc) = centroids(cc,1);
            centroidsByFrame_Y(ii,cc) = centroids(cc,2);
            % also get the frame ready for display, if needed
            if p.Results.displayMode && ~p.Results.useParallel
                thisFrame = insertShape(thisFrame,'FilledCircle', [centroids(cc,1),centroids(cc,2),2], 'Color','red');
            end
        end  
    end
    
    % display the frame if requested
    if p.Results.displayMode && ~p.Results.useParallel
        imshow(thisFrame,'Border', 'tight', 'InitialMagnification', 200)
    end
end

%% Refine data
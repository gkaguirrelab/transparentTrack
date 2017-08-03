function [glintData] = trackGlint(grayVideoName, glintFileName, varargin)
% function [glintData] = trackGlint(grayVideoName, glintFileName, varargin)
%
% This function tracks the glint using the circle patch + direct ellipse
% fitting approach.
%
% The routine first identifies the likely location of the pupil, and then
% dilates this to create a search region. The glint is then found within.
%
% Several parameters control the glint search. The routine is robust and
% rarely should adjustment in the parameters be needed to process a video
% successfully.
%
% Output
%	glintData - structure with fields that contain the X and Y location of
%       the center of the glint (in units of pixels), as well as the
%       ellipseFittingError that is set to true if the X and Y position of
%       the glint is based on the circle patch only.
%
%
% Input (required)
%	grayVideoName - full path to the video in which to track the glint.
%   	A grayscale video is expected.
%   glintFileName - full path to the matFile in which to save the glint
%       results.
%
% Options (analysis)
%	glintGammaCorrection - gamma correction to be applied in current
%       frame. An extremely high value will make almost all the frame black
%       and only big bright spots will be white. This reduces the
%       possibility of confusing the glint with some other smaller bright
%       spot (default 5.5, decrease if no glint is found)
%   glintCircleThresh - relative threshold value to locate the glint for
%       circle fitting. The high number used reflects the fact that the
%       glint should be the brightest point within the search region.
%	glintRange - radius range for circle fitting of the glint (in
%       pixels)
%   glintEllipseThresh - threshold value to locate the glint for
%     	ellipse fitting (default 0.9)
%   pupilRange - pupil range initialization for more accurate glint
%       tracking in the circle patch step
%   pupilCircleThresh - pupil threshold initialization for more accurate
%       glint tracking in the circle patch step
%
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
p.addParameter('glintGammaCorrection', 5.5, @isnumeric);
p.addParameter('glintCircleThresh', 0.999, @isnumeric);
p.addParameter('glintRange', [10 30], @isnumeric);
p.addParameter('glintEllipseThresh', 0.9, @isnumeric);
p.addParameter('pupilRange', [30 90], @isnumeric);
p.addParameter('pupilCircleThresh', 0.06, @isnumeric);
p.addParameter('frameMask',[] , @isnumeric);

% Optional display params
p.addParameter('verbosity','none',@ischar);
p.addParameter('displayMode',false,@islogical);

% CircleFit routine params. Defined here for transparency
p.addParameter('pupilOnly', false, @islogical);
p.addParameter('glintOut', 0.15, @isnumeric);
p.addParameter('dilateGlint', 6, @isnumeric);
p.addParameter('imfindcirclesSensitivity', 0.99, @isnumeric);
p.addParameter('rangeAdjust', 0.05, @isnumeric);

% Optional flow control params
p.addParameter('nFrames',[],@isnumeric);
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
% read the video into memory, adjusting gamma if needed
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


%% Track the glint

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

% Initialize glint variables
glintData_X = nan(nFrames,1);
glintData_Y = nan(nFrames,1);
glintData_ellipseFittingError = nan(nFrames,1);

%loop through frames
parfor (ii = 1:nFrames, nWorkers)
    
    % increment the progress bar
    if strcmp(p.Results.verbosity,'full') && mod(ii,round(nFrames/50))==0
        fprintf('\b.\n');
    end
    
    % explicitly define these variables so that parfor does not worry that
    % they are expected to carry across loops
    Egi=0;
    Eg=0;
    
    % get the frame
    thisFrame = squeeze(grayVideo(:,:,ii));
    
    % apply a frame mask if requested
    if ~isempty (p.Results.frameMask)
        thisFrame((1:p.Results.frameMask(1)),:) = 220;
        thisFrame((end - p.Results.frameMask(1):end),:) = 220;
        thisFrame(:, (1:p.Results.frameMask(2))) = 220;
        thisFrame(:, (end - p.Results.frameMask(2):end)) = 220;
    end
    
    % track with circles (using the default options)
    [~,~,~, gCenters, gRadii,~, ~, ~] = findGlintAndPupilCircles(thisFrame, ...
        p.Results.pupilCircleThresh, ...
        p.Results.glintCircleThresh, ...
        p.Results.pupilRange, ...
        p.Results.glintRange,...
        p.Results.pupilOnly,p.Results.glintOut,p.Results.dilateGlint,p.Results.imfindcirclesSensitivity,p.Results.rangeAdjust);
    
    % if a glint was present, refine the location
    if ~isempty(gCenters)
        % getGlintPerimeter
        [binG] = getGlintPerimeter (thisFrame, gCenters, gRadii, p.Results.glintEllipseThresh);
        % Fit ellipse to glint
        [Yg, Xg] = ind2sub(size(binG),find(binG));
        % turn of warnings for singular matrix
        origWarnState = warning();
        warning('off','MATLAB:singularMatrix');
        warning('off','MATLAB:illConditionedMatrix');
        try
            Eg = ellipse_im2ex(ellipsefit_direct(Xg,Yg));
            if isreal(Eg)
                glintData_X(ii) = Eg(1);
                glintData_Y(ii) = Eg(2);
            else
                glintData_X(ii)= gCenters(1,1);
                glintData_Y(ii) = gCenters(1,2);
                glintData_ellipseFittingError(ii) = 1;
            end
        catch % "Index exceeds matrix dimensions" for ellipsefit_direct
            glintData_X(ii)= gCenters(1,1);
            glintData_Y(ii) = gCenters(1,2);
            glintData_ellipseFittingError(ii) = 1;
        end
        warning(origWarnState);
    end
    
    % display
    if p.Results.displayMode && ~p.Results.useParallel
        if ~isnan(glintData_X(ii)) && glintData_ellipseFittingError(ii)==1
            dispFrame = insertShape(thisFrame,'FilledCircle', [glintData_X(ii),glintData_Y(ii), 2],'Color','yellow');
        elseif ~isnan(glintData_X(ii)) && isnan(glintData_ellipseFittingError(ii))
            dispFrame = insertShape(thisFrame,'FilledCircle', [glintData_X(ii),glintData_Y(ii), 2],'Color','red');
        else
            dispFrame = thisFrame;
        end
        
        imshow(dispFrame,'Border', 'tight')
    end
    
end

%% Clean up, save, and close

% gather loop variables into a structure
glintData.X = glintData_X;
glintData.Y = glintData_Y;
glintData.ellipseFittingError = glintData_ellipseFittingError;

% add a meta field with analysis details
glintData.meta = p.Results;

% save out a mat file with the glint tracking data
if ~p.Results.displayMode
    save (glintFileName, 'glintData')
else
    close(figureHandle);
end

% report completion of analysis
if strcmp(p.Results.verbosity,'full')
    toc
    fprintf('\n');
end

% Delete the parallel pool
if strcmp(p.Results.verbosity,'full')
    tic
    fprintf(['Closing parallel pool. Started ' char(datetime('now')) '\n']);
end
if p.Results.useParallel
    poolObj = gcp;
    if ~isempty(poolObj)
        delete(poolObj);
    end
end
if strcmp(p.Results.verbosity,'full')
    toc
    fprintf('\n');
end

end % main function


function [binG] = getGlintPerimeter (I, gCenters, gRadii, glintEllipseThresh)
% create a mask from circle fitting parameters (note: glint
% is already dilated
glintMask = zeros(size(I));
glintMask = insertShape(glintMask,'FilledCircle',[gCenters(1,1) gCenters(1,2) gRadii(1)],'Color','white');
glintMask = im2bw(glintMask);

% apply mask to grey image
maskedGlint = immultiply(I,glintMask);

% convert back to gray
gI = uint8(maskedGlint);

% Binarize glint
binG  = ones(size(gI));
binG(gI<quantile(double(I(:)),glintEllipseThresh)) = 0;

% get perimeter of glint
binG = bwperim(binG);

end % getGlintPerimeter

function [perimeter] = extractPupilPerimeter(grayVideoName, perimeterFileName, varargin)
% function [perimeterParams] = extractPupilPerimeter(grayVideoName, perimeterVideoName,varargin)
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
%	pupilEllipseThresh - threshold value to locate the glint for
%       ellipse fitting (default 0.9, typical range [0.8 0.98])
%   glintCircleThresh - DEFINE HERE
%   glintRange - DEFINE HERE
%   maskBox - DEFINE HERE
%   smallObjThresh - DEFINE HERE
% 
% Options (display)
%   verbosity - controls console status updates
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


%% Parse input and define variables
p = inputParser; p.KeepUnmatched = true;

% required input
p.addRequired('grayVideoName',@isstr);
p.addRequired('perimeterFileName',@isstr);

% Optional analysis params
p.addParameter('gammaCorrection', 1, @isnumeric);
p.addParameter('pupilCircleThresh', 0.06, @isnumeric);
p.addParameter('pupilRange', [30 90], @isnumeric);
p.addParameter('pupilEllipseThresh', 0.945, @isnumeric);
p.addParameter('glintCircleThresh', 0.999, @isnumeric);
p.addParameter('glintRange', [10 30], @isnumeric);
p.addParameter('maskBox', [4 30], @isnumeric);
p.addParameter('smallObjThresh', 500, @isnumeric);

% Optional display params
p.addParameter('verbosity','none',@ischar);

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
% read the video into memory, adjusting gamma if needed
for ii = 1:nFrames
    thisFrame = readFrame(videoInObj);
    thisFrame = imadjust(thisFrame,[],[],p.Results.gammaCorrection);
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


%% Extract pupil perimeter

% alert the user
if strcmp(p.Results.verbosity,'full')
    tic
    fprintf(['Extracting pupil perimeter. Started ' char(datetime('now')) '\n']);
    fprintf('| 0                      50                   100%% |\n');
    fprintf('.\n');
end

% initialize variable to hold the perimeter data
perimeter_data = zeros(videoSizeY,videoSizeX,nFrames,'uint8');

% structuring element for pupil mask size
sep = strel('rectangle',p.Results.maskBox);

% loop through gray frames
parfor (ii = 1:nFrames, nWorkers)
    % increment the progress bar
    if strcmp(p.Results.verbosity,'full') && mod(ii,round(nFrames/50))==0
        fprintf('\b.\n');
    end
    
    % get the frame
    thisFrame = squeeze(grayVideo(:,:,ii));
    
    % perform an initial search for the pupil with circleFit 
    [pCenters, pRadii,~,~,~,~, pupilRange, ~] = ...
        circleFit(thisFrame,...
        p.Results.pupilCircleThresh,...
        p.Results.glintCircleThresh,...
        p.Results.pupilRange,...
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
        perimeter_data(:,:,ii) = thisFrame;
    else
        thisFrame = im2uint8(zeros(size(thisFrame)));
        perimeter_data(:,:,ii) = thisFrame;
    end
            
end % loop through gray frames

%% Clean up and close

% save mat file with the video and analysis details
perimeter.data = perimeter_data;
perimeter.meta = p.Results;
save(perimeterFileName,'perimeter');

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

end % function

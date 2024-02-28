function perimeter = findPupilPerimeter(grayVideoName, perimeterFileName, varargin)
% Threshold video frames to find the pupil perimeter
%
% Syntax:
%  perimeter = findPupilPerimeter(grayVideoName, perimeterFileName)
%
% Description:
%	An initial search for the pupil border is performed with the local
%   function 'findPupilCircle'. If a candidate circle is found, the region
%   is dilated. We then binarize the resulting "patch" image with a user
%   set threshold, and extract the perimeter of the bigger region surviving
%   the thresholding process (believed to be the pupil).
%
% Inputs:
%	grayVideoName -       - Full path to the gray video to track
%	perimeterFileName     - Full path to the output .mat file
%
% Optional key/value pairs (display and I/O):
%  'verbose'              - Logical. Default false.
%  'displayMode'          - When set to true, displays the results of the
%                           boundary extraction and does not save a video
%
% Optional key/value pairs (flow control)
%  'nFrames'              - Analyze fewer than the total number of frames
%  'startFrame'           - Which frame to start on
%
% Optional key/value pairs (environment)
%  'tbSnapshot'           - This should contain the output of the
%                           tbDeploymentSnapshot performed upon the result
%                           of the tbUse command. This documents the state
%                           of the system at the time of analysis.
%  'timestamp'            - AUTOMATIC; The current time and date
%  'username'             - AUTOMATIC; The user
%  'hostname'             - AUTOMATIC; The host
%
% Optional key/value pairs (analysis)
%  'pupilGammaCorrection' - Gamma correction to be applied to the video
%                           frames (default 1, typical values between 0.5
%                           and 1.8)
%  'pupilFrameMask'       - This option will add a mask on the original
%                           gray video, framing it by [nRows nColumns] on
%                           the borders symmetrically or by [nRowsTop
%                           nColumnsRight nRowsBottom nColumnsLeft]. This
%                           is particularly useful for size calibration
%                           videos in which appear partial black dots that
%                           may throw off the circle finding mechanism.
%                           fitting (default 0.06, typical range [0.04
%                           0.09])
%  'maskBox'              - This is the proportion to dilate the pupil
%                           masked region in the vertical and horizontal
%                           directions respectively. A value of 1 will
%                           result in no dilation in that direction. A
%                           value of 2 will result in a masked region that
%                           is twice the size of the pupil radius.
%  'frameMaskValue'       - The image value that is assigned to the region
%                           that is masked by frameMask. This should be a
%                           gray that is neither pupil nor glint.
%  'smallObjThresh'       - Maximum size of small objects to be removed to
%                           clean up the pupil perimeter.
%  'expandPupilRange'     - This option controls what the routine does when
%                           it does not find a circle within the initial
%                           range. If true, the routine will expand the
%                           pupilRange and look again. If false and no
%                           circle is found initially, the routine will
%                           move onto the next frame. Consider setting this
%                           to false when the routine finds circles that
%                           don't correspond to the pupil.
%
% Optional key/value pairs (used in the local function findPupilCircle)
%  'pupilCircleThresh'    - The threshold used to binarize the image.
%                           Ranges from 0 to 1, with 0 corresponding to
%                           black.
%  'adaptivePupilRangeFlag' - Logical. If set to true, then the pupil range
%                           is updated on each frame based upon the prior
%                           frame.
%  'pupilRange'           - Initial radius range in pixels for circle
%                           fitting of the pupil (default [20 100]). This
%                           value may be dynamically updated.
%  'imfindcirclesSensitivity' - A parameter (ranging from 0-1) that is used
%                           in the call to the matlab function
%                           imfindcircles as part of the key-value pair
%                           'Sensitivity'. We set the default value to
%                           0.99, meaning that the function will be very
%                           sensitive to circles, even those that are
%                           partially obscured.
%  'rangeAdjust'          - For each frame, a circle is fit to the pupil.
%                           The size of the found pupil is used to update
%                           the pupilRange values for the next frame. The
%                           rangeAdjust parameter defines the proportion
%                           above and below this target size that is used
%                           to search on the next frame
%  'pickLargestCircle'    - This routine uses imfindcircles to identify the
%                           likely pupil. However, it can return more than
%                           1 circle that it finds within the
%                           pupilMaskFrame. By default it sorts according
%                           to most circular, but sometimes the pupil is
%                           more likely to be the biggest circle it finds.
%                           The default is set to false.
%
% Outputs:
%   perimeter             - Structure with a 'data' field that contains the
%                           subfields Xp and Yp, which contain the sub
%                           indices of the points the constitute the pupil
%                           perimeter. There is a 'meta' field with
%                           analysis and environment params
%


%% Parse input and define variables
p = inputParser; p.KeepUnmatched = true; p.PartialMatching = false;

% Required input
p.addRequired('grayVideoName',@isstr);
p.addRequired('perimeterFileName',@isstr);

% Optional display and I/O params
p.addParameter('verbose',false,@islogical);
p.addParameter('displayMode',false,@islogical);

% Optional flow control params
p.addParameter('nFrames',Inf,@isnumeric);
p.addParameter('startFrame',1,@isnumeric);

% Optional environment parameters
p.addParameter('tbSnapshot',[],@(x)(isempty(x) | isstruct(x)));
p.addParameter('timestamp',char(datetime('now')),@ischar);
p.addParameter('username',char(java.lang.System.getProperty('user.name')),@ischar);
p.addParameter('hostname',char(java.net.InetAddress.getLocalHost.getHostName),@ischar);

% Optional analysis params
p.addParameter('adaptivePupilRangeFlag', true, @islogical);
p.addParameter('pupilGammaCorrection', 0.75, @isnumeric);
p.addParameter('pupilFrameMask', [], @isnumeric);
p.addParameter('maskBox', [2 2], @isnumeric);
p.addParameter('nOtsu', 2, @isnumeric);
p.addParameter('frameMaskValue', 220, @isnumeric);
p.addParameter('smallObjThresh', 400, @isnumeric);
p.addParameter('expandPupilRange', true, @islogical);
p.addParameter('pickLargestCircle', false, @islogical);
p.addParameter('glintFileName', [], @ischar);
p.addParameter('recentGlintRange', 5, @isnumeric);


% Optional findPupilCircle routine params. Defined here for transparency
p.addParameter('pupilCircleThresh', 0.06, @isnumeric);
p.addParameter('pupilRange', [20 100], @isnumeric);
p.addParameter('imfindcirclesSensitivity', 0.99, @isnumeric);
p.addParameter('rangeAdjust', 0.05, @isnumeric);

% parse
p.parse(grayVideoName, perimeterFileName, varargin{:})


% Prepare the video object
videoInObj = videoIOWrapper(grayVideoName,'ioAction','read');

% get number of frames
if p.Results.nFrames == Inf
    if isprop(videoInObj,'NumFrames')
        nFrames = videoInObj.NumFrames;
    else
        nFrames = floor(videoInObj.Duration*videoInObj.FrameRate);
    end
else
    nFrames = p.Results.nFrames;
end


% get video dimensions
videoSizeX = videoInObj.Width;
videoSizeY = videoInObj.Height;

% Get the number of regions into which the image should be segmented
nOtsu = p.Results.nOtsu;

% initialize variable to hold the perimeter data
grayVideo = zeros(videoSizeY,videoSizeX,nFrames,'uint8');
cc = 0;
for ii = p.Results.startFrame:p.Results.startFrame+nFrames-1
    cc = cc+1;
    thisFrame = read(videoInObj,ii);
    thisFrame = imadjust(thisFrame,[],[],p.Results.pupilGammaCorrection);
    grayVideo(:,:,cc) = rgb2gray (thisFrame);
end
% end
% close the video object
% figure out which frames we're looping through
startFrame = p.Results.startFrame;
if p.Results.nFrames == Inf
    endFrame = floor(videoInObj.Duration*videoInObj.FrameRate);
else
    endFrame = startFrame + p.Results.nFrames - 1;
end
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
if p.Results.verbose
    tic
    fprintf(['Extracting pupil perimeter. Started ' char(datetime('now')) '\n']);
    fprintf('| 0                      50                   100%% |\n');
    fprintf('.');
end

% initialize variable to hold the perimeter data
perimeter = struct();
perimeter.size = [size(thisFrame,1) size(thisFrame,2)];
perimeter.data = cell(nFrames,1);

% Initialize the pupilRange with the parameter value. This value is updated
% as we progress through the frames
pupilRange = p.Results.pupilRange;


% loop through gray frames

for ii = 1:(endFrame-startFrame+1)
    
    if p.Results.displayMode && strcmp(get(figureHandle,'currentchar'),' ')
        close(figureHandle)
        return
    end
    
    % increment the progress bar
    if p.Results.verbose && mod(ii-p.Results.startFrame+1,round(nFrames/50))==0
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

    % If we are not using an adaptive pupil range search, just
    % re-initialize with the passed parameter
    if ~p.Results.adaptivePupilRangeFlag
        pupilRange = p.Results.pupilRange;
    end
        
    % if provided with glintFileName, only look for perimeter if glints are in
    % the present frame. If no glints are present, this frame is a blink and
    % any frame would be a bad frame
    if ~isempty(p.Results.glintFileName)
        load(p.Results.glintFileName);
        frameNumber = startFrame + ii - 1;
        % if this frame is a blink, don't bother looking for a circle
        if sum(isnan(glintData.X(frameNumber,:))) == 2
            performFindPupilCircle = false;
            pCenters = [];
        else
            performFindPupilCircle = true;
        end
        % if the previous frame was a blink, start search from initial
        % pupil range (rather than potentially spurious range that would be
        % problematic
        if frameNumber - p.Results.recentGlintRange < 1
            recentGlints = glintData.X(1:frameNumber,:);
        else
            recentGlints = glintData.X(frameNumber - p.Results.recentGlintRange:frameNumber, :);
        end
        
        % if any of the recent frames were blinks, restore pupilRange to the
        % initial pupil range
        if sum(sum(isnan(recentGlints)))>0
            pupilRange = p.Results.pupilRange;
        end
        
        
    else
        performFindPupilCircle = true;
    end
    
    if performFindPupilCircle
        
        % perform an initial search for the pupil with findGlintAndPupilCircles. Also extract
        % glint location and size information for later use.
        [pCenters, pRadii,~,pupilRange] = ...
            findPupilCircle(thisFrame,...
            p.Results.pupilCircleThresh,...
            pupilRange,...
            p.Results.imfindcirclesSensitivity,p.Results.rangeAdjust, 'displayMode', p.Results.displayMode);
        
        if (p.Results.expandPupilRange)
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
                    p.Results.imfindcirclesSensitivity,p.Results.rangeAdjust, 'displayMode', p.Results.displayMode);
                if isempty(pCenters) % still no circle? Try 100% increase
                    candidateRange = [ceil(initialPupilRange(1)/2) round(initialPupilRange(2)*2)];
                    if candidateRange(1) < p.Results.pupilRange(1) || candidateRange(2) > p.Results.pupilRange(2)
                        candidateRange = p.Results.pupilRange;
                    end
                    [pCenters, pRadii,~,pupilRange] = ...
                        findPupilCircle(thisFrame,...
                        p.Results.pupilCircleThresh,...
                        candidateRange,...
                        p.Results.imfindcirclesSensitivity,p.Results.rangeAdjust, 'displayMode', p.Results.displayMode);
                    if isempty(pCenters) % STILL no circle? Give up and restore initialPupilRange
                        pupilRange = initialPupilRange;
                    end
                end
            end
        end
    end
    
    % If a pupil circle patch was ultimately found, get the perimeter, else
    % write out a zero-filled frame
    if ~isempty(pCenters)
        
        if p.Results.pickLargestCircle
            [ pRadii, sortIndices ] = sort(pRadii, 'descend');
            pCenters = pCenters(sortIndices,:);
            
        end
        
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
        otsuThresh = multithresh(maskedPupilNaN,nOtsu);
        
        % restore warning state
        warning(warningState);
        
        % Set the glint and iris to zero, and the pupil to unity
        vec = zeros(1,nOtsu+1);
        vec(end)=1;
        binP = imquantize(maskedPupil, otsuThresh, vec);
        
        % remove small objects
        binP = bwareaopen(binP, p.Results.smallObjThresh);
        
        % fill the holes
        binP = imfill(binP,'holes');
        
        % get perimeter of object
        binP = bwperim(binP);
        
        % save the perimeter
        [perimeter.data{ii}.Yp, perimeter.data{ii}.Xp] = ind2sub(size(binP),find(binP));
    else
        perimeter.data{ii}.Yp = [];
        perimeter.data{ii}.Xp = [];
    end
    
    if p.Results.displayMode
        displayFrame=thisFrame;
        if ~isempty(perimeter.data{ii}.Xp)
            displayFrame(sub2ind(size(thisFrame),perimeter.data{ii}.Yp,perimeter.data{ii}.Xp))=255;
        end
        axes('units','norm','outerposition',[0.5 0 0.5 1],'position',[0.5 0 0.5 1]) %right axes
        imshow(displayFrame, 'Border', 'tight');
        if isempty(perimeter.data{1}.Xp)
            string = 'No pupil found';
        else
            string = (['Frame ', num2str(p.Results.startFrame)]);
        end
        text(350, 200, string);
        
        set(gcf, 'Position', [274 462 1448 430]);
    end
    
end % loop through gray frames

%% Clean up and close

% save mat file with the video and analysis details
perimeter.meta = p.Results;
if ~p.Results.displayMode
    save(perimeterFileName,'perimeter','-v7.3');
else
    %close(figureHandle);
end

% report completion of analysis
if p.Results.verbose
    fprintf('\n');
    toc
    fprintf('\n');
end


end %  main function



%% LOCAL FUNCTIONS

function [pCenters, pRadii,pMetric, pupilRange] = findPupilCircle(I,pupilCircleThresh,pupilRange,imfindcirclesSensitivity,rangeAdjust, varargin)
% This function is used for pupil circle fitting.

%% parse input and define variables
p = inputParser;
% required input
p.addRequired('I');
p.addRequired('pupilCircleThresh', @isnumeric)
p.addRequired('pupilRange',@isnumeric);
p.addRequired('imfindcirclesSensitivity', @isnumeric);
p.addRequired('rangeAdjust', @isnumeric);

% optional input
p.addParameter('displayMode',false,@islogical);


% parse
p.parse(I,pupilCircleThresh,pupilRange,imfindcirclesSensitivity,rangeAdjust, varargin{:});

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

if p.Results.displayMode
    axes('units','norm','outerposition',[0 0 0.5 1],'position',[0 0 0.5 1]) %left axes
    imshow(binP)
    if ~isempty(pCenters)
        hold on; viscircles(pCenters, pRadii);
        viscircles(pCenters(1,:), pRadii(1), 'Color', 'b');
        text(350, 200, [' Radius = ', num2str(pRadii(1))]);
    end
end

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

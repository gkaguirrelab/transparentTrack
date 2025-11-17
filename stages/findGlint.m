function [glintData] = findGlint(grayVideoName, glintFileName, varargin)
% Identifies one or more glints in the frames of an IR video of the eye
%
% Syntax:
%  [glintData] = findGlint(grayVideoName, glintFileName)
%
% Description:
%   This function tracks one or more glints in an IR video using a simple
%   thresholding and region property identification approach.
%
%   Every frame is first corrected with a gamma > 1, so that glints and
%   other large bright spot are enhanced. The image is then binarized with
%   a relatively high threshold, so that only the bright spots and their
%   immediate surroundings remain as non-zero values. The centroids of each
%   surviving region in the binary image are extracted using the MATLAB
%   function "regionprops". The centroid location is weighted with the
%   actual brightness value of each pixel in the gray gamma-corrected
%   image.
%
%   After all centroid locations are extracted, data is refined according
%   to the expected number of glints and average centroid location
%   throughout the video. At the moment, the routine is able to process
%   videos in which either 1 or 2 glints need to be tracked.
%
%   In the single glint case, we calculate the median location of the glint
%   from those frames that return as just a single centroid. For the
%   frames in which more than the expected number of glints is found, we
%   use the median value of the "good centroids" location to assess which
%   of the bright regions identified is indeed the desired glint.
%
%   In the double glint case, we assume that the couple of glints has a
%   "main direction" that the user needs to declare. That is: if the 2
%   light sources that produce the glint are vertically spaced in the real
%   world, the main direction will be 'y' on the video. Conversely, if they
%   are horizontally spaced, the main direction will be 'x'. As in the
%   single glint case, we start by locating the frames that present the
%   desired number of bright spot (2).Assuming that in most of those frames
%   we are tracking the correct glints, we compute a positive "glint
%   vector" using the tracked centroids coordinate. That is to derive both
%   an orientation (that is used to sort the glint couples consistently
%   from one frame to the other) and a median glint vector length. In the
%   subsequent step, the median glint vector lenght is used to narrow down
%   the most likely candidates glint in frames where more than the desired
%   number of centroids was tracked.
%
%   Both for the single and double glint case, in frames where less than
%   the desired number of glints is located, we assume that no reliable
%   glint information was available (e.g. subject was blinking) and the
%   corresponding glint information will be set as NaNs.
%
% Notes:
%   Coordinate system - the function "regionprops" will save the centroids
%   in world coordinates (origin top left corner of the frame, xlim = [0
%   horizontalRes], ylim = [0 verticalRes], therefore, the glint data will
%   also be expressed in world coordinates.
%
% Inputs:
%	grayVideoName         - Full path to the video in which to track the
%                           glint. A grayscale video is expected
%   glintFileName         - Full path to the glint file
%
% Optional key/value pairs (display and I/O):
%  'verbose'              - Logical. Default false.
%  'displayMode'          - If set to true, a continuously updated video
%                           displays the glint fitting. This is slow but
%                           may be useful while setting analysis params.
%
% Optional key/value pairs (flow control)
%  'nFrames'              - Analyze fewer than the total number of frames.
%  'startFrame'           - First frame from which to start the analysis.
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
%  'numberOfGlints'       - Desired number of glints to find
%  'glintsMainDirection'  - In the two glints case, this is the main
%                           direction in which the two glints are
%                           consistently spaced.
%  'glintGammaCorrection' - Gamma correction to be applied in current
%                           frame. An extremely high value will make almost
%                           all the frame black and only big bright spots
%                           will be white. This reduces the possibility of
%                           confusing the glint with some other smaller
%                           bright spot (default 1.5, decrease if no glint
%                           is found, increase if too many glints are
%                           found)
%  'glintThreshold' 	  - Threshold value to binarize the glint gray
%                           image. Should be set to preserve both the
%                           glints and any "halo" around them.(default
%                           value 0.8)
%  'glintFrameMask'       - Add a mask on the original gray video, framing
%                           it by [nRows nColumns] on the borders
%                           symmetrically or by [nRowsTop nColumnsRight
%                           nRowsBottom nColumnsLeft].
%  'frameMaskValue'       - The image value that is assigned to the region
%                           that is masked by frameMask. This should be a
%                           gray that is neither pupil nor glint.
%  'centroidsAllocation'  - Max number of centroids to be saved in memory
%  'medianGlintVector'    - Optionally specify the true x and y distances
%                           between the glints. If unspecified, the routine
%                           will determine these values on its own.
%
% Output
%	glintData             - Structure with fields that contain the X and Y
%                           location of the center of the glint (in units
%                           of pixels), and a meta field with additional
%                           analysis params and intermediate results.
%
% Examples:
%{
    videoInFileName = ('/Users/aguirre/Aguirre-Brainard Lab Dropbox/Geoffrey Aguirre/FLIC_analysis/lightLogger/scriptedIndoorOutdoor/FLIC_2004/gazeCalibration/temporalFrequency/FLIC_2004_gazeCalibration_tf_session1/P.avi');
    perimeterFileName = ('/Users/aguirre/Aguirre-Brainard Lab Dropbox/Geoffrey Aguirre/FLIC_analysis/lightLogger/scriptedIndoorOutdoor/FLIC_2004/gazeCalibration/temporalFrequency/FLIC_2004_gazeCal_perimeter.mat');
    glintFileName = ('/Users/aguirre/Aguirre-Brainard Lab Dropbox/Geoffrey Aguirre/FLIC_analysis/lightLogger/scriptedIndoorOutdoor/FLIC_2004/gazeCalibration/temporalFrequency/FLIC_2004_gazeCal_glint.mat');
    findGlint(videoInFileName, glintFileName,'perimeterFileName',perimeterFileName,'verbose',true,'startFrame',12000,'glintFrameMask',[125,150,200,150]);
%}

%% parse input and define variables
p = inputParser; p.KeepUnmatched = true; p.PartialMatching = false;

% Required
p.addRequired('grayVideoName',@isstr);
p.addRequired('glintFileName',@isstr);

% Optional display and I/O params
p.addParameter('verbose',false,@islogical);
p.addParameter('displayMode',false,@islogical);

% Optional flow control params
p.addParameter('nFrames',Inf,@isnumeric);
p.addParameter('startFrame',1,@isnumeric);

% Optional environment params
p.addParameter('tbSnapshot',[],@(x)(isempty(x) | isstruct(x)));
p.addParameter('timestamp',char(datetime('now')),@ischar);
p.addParameter('username',char(java.lang.System.getProperty('user.name')),@ischar);
p.addParameter('hostname',char(java.net.InetAddress.getLocalHost.getHostName),@ischar);

% Optional analysis params
p.addParameter('numberOfGlints', 1, @isnumeric);
p.addParameter('glintsMainDirection', 'y',@ischar);
p.addParameter('glintGammaCorrection', 1.5, @isnumeric);
p.addParameter('glintThreshold', 0.8, @isnumeric);
p.addParameter('glintFrameMask',[] , @isnumeric);
p.addParameter('frameMaskValue', 30, @isnumeric);
p.addParameter('centroidsAllocation', 5, @isnumeric);
p.addParameter('threshold',[], @isnumeric);
p.addParameter('medianGlintVector',[],@isnumeric);
p.addParameter('perimeterFileName','',@ischar);
p.addParameter('removeIsolatedGlints',false,@islogical);



% parse
p.parse(grayVideoName, glintFileName, varargin{:})


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

%% Initialize variables
% Initialize glint variables
glintData_X = nan(nFrames,p.Results.numberOfGlints);
glintData_Y = nan(nFrames,p.Results.numberOfGlints);

% initialize centroid variable according to the centroidsAllocation value
centroidsByFrame_X = nan(nFrames,p.Results.centroidsAllocation);
centroidsByFrame_Y = nan(nFrames,p.Results.centroidsAllocation);


%% Find all centroids

% Detect display mode
if p.Results.displayMode
    fprintf('** DISPLAY MODE **\n')
    fprintf('Results will not be saved.\n')
    
    % create a figure for display
    figureHandle=figure();
    
    % we will monitor the currentchar for a 'q'
    set(figureHandle,'currentchar','?')
end

% alert the user
if p.Results.verbose
    tic
    fprintf(['Tracking the glint. Started ' char(datetime('now')) '\n']);
    fprintf('| 0                      50                   100%% |\n');
    fprintf('.\n');
end

% loop through frames
for ii = 1:nFrames
    if p.Results.displayMode && strcmp(get(figureHandle,'currentchar'),' ')
        close(figureHandle)
        return
    end
    
    % increment the progress bar
    if p.Results.verbose && mod(ii,round(nFrames/50))==0
        fprintf('\b.\n');
    end

    if ii<p.Results.startFrame
        readFrame(videoInObj);
        while videoInObj.last_frame_read<ii
        end
        continue
    end
    
    % get the frame
    thisFrame = readFrame(videoInObj);
    thisFrame = imadjust(thisFrame,[],[],p.Results.glintGammaCorrection);
    thisFrame = rgb2gray (thisFrame);
    
    % apply a frame mask if required
    if ~isempty (p.Results.glintFrameMask)
        if length(p.Results.glintFrameMask) == 2
            thisFrame((1:p.Results.glintFrameMask(1)),:) = p.Results.frameMaskValue;
            thisFrame((end - p.Results.glintFrameMask(1):end),:) = p.Results.frameMaskValue;
            thisFrame(:, (1:p.Results.glintFrameMask(2))) = p.Results.frameMaskValue;
            thisFrame(:, (end - p.Results.glintFrameMask(2):end)) = p.Results.frameMaskValue;
        elseif length(p.Results.glintFrameMask) == 4
            thisFrame((1:p.Results.glintFrameMask(1)),:) = p.Results.frameMaskValue; %top
            thisFrame(:, (end - p.Results.glintFrameMask(2):end)) = p.Results.frameMaskValue; %left
            thisFrame((end - p.Results.glintFrameMask(3):end),:) = p.Results.frameMaskValue; %bottom
            thisFrame(:, (1:p.Results.glintFrameMask(4))) = p.Results.frameMaskValue; %right
        else
            error ('invalid frameMask parameter. Frame mask must be defined as [nRows nColumns] or as [nRowsTop nColumnsRight nRowsBottom nColumnsLeft]')
        end
    end
    
    % binarize glint image according to glintThreshold
    binG = imbinarize(thisFrame, p.Results.glintThreshold);
    
    % get the weighted centroids for all the surviving bright spots in the
    % image. The "weight" is given by the actual gray value in the frame,
    % so that the brightest glint pixels are more relevant to the overall
    % glint position.
    stats = regionprops(binG, thisFrame, 'WeightedCentroid', 'MajorAxisLength','MinorAxisLength');
    
    % order glints by size
    centroidSizes = [stats.MajorAxisLength];
    [~, sizeOrder] = sort(centroidSizes, 'descend');
    
    % if centroids were found in this frame, save them out.
    if ~isempty(stats)
        centroidCounter = 1;
        for jj = sizeOrder
            centroids(centroidCounter,:) = stats(jj).WeightedCentroid;
            centroidCounter = centroidCounter + 1;
        end
        clear stats
        for cc = 1: min(size(centroids,1),p.Results.centroidsAllocation)
            centroidsByFrame_X(ii,cc) = centroids(cc,1);
            centroidsByFrame_Y(ii,cc) = centroids(cc,2);
            % also get the frame ready for display, if needed
            if p.Results.displayMode
                thisFrame = insertShape(thisFrame,'FilledCircle', [centroids(cc,1),centroids(cc,2),2], 'Color','red');
            end
        end
        clear centroids
    end
    
    % display the frame if requested
    if p.Results.displayMode
        imshow(thisFrame,'Border', 'tight', 'InitialMagnification', 200)
    end

end


%% Get glint data

% get the number of centroids found for each frame
centroidsInEachFrame = sum(~isnan(centroidsByFrame_X),2);

% first, save out data for the frames with the expected amount of glints
framesWithExpectedCentroids = find (centroidsInEachFrame==p.Results.numberOfGlints);

switch p.Results.numberOfGlints
    case 1 % this case is simple and does not require clustering of the centroids
        for ii = 1: length(framesWithExpectedCentroids)
            glintData_X(framesWithExpectedCentroids(ii)) = centroidsByFrame_X(framesWithExpectedCentroids(ii),1);
            glintData_Y(framesWithExpectedCentroids(ii)) = centroidsByFrame_Y(framesWithExpectedCentroids(ii),1);
        end
        
    case 2
        % generate an array with a glint-to-glint vector from the frames
        % with expected glints.
        if isempty(p.Results.medianGlintVector)
            if ~isempty(framesWithExpectedCentroids)
                for ii = 1: length(framesWithExpectedCentroids)
                    unsortedGlintVector(ii,1) = centroidsByFrame_X(framesWithExpectedCentroids(ii),1) - centroidsByFrame_X(framesWithExpectedCentroids(ii),2);
                    unsortedGlintVector(ii,2) = centroidsByFrame_Y(framesWithExpectedCentroids(ii),1) - centroidsByFrame_Y(framesWithExpectedCentroids(ii),2);
                end
                
                % get median length components of the glint vector
                medianGlintVector = [nanmedian(abs(unsortedGlintVector(:,1))) nanmedian(abs(unsortedGlintVector(:,2)))];
            else % in case all frames do not have the expected number of glints
                unsortedGlintVector = [];
                % we're going to figure out all of the possible distances
                % between glints
                xLengths= [];
                yLengths = [];
                for ii = 1:length(centroidsByFrame_X)
                    newXLengths = pdist(centroidsByFrame_X(ii,:)');
                    xLengths = [xLengths, newXLengths];
                    
                    newYLengths = pdist(centroidsByFrame_Y(ii,:)');
                    yLengths = [yLengths, newYLengths];
                end
                
                % then figure out which is most common, and use that as our
                % medianGlintVector
                [f, x] = ksdensity(xLengths);
                [~, maxIndex] = max(f);
                medianGlintVector(1) = x(maxIndex);
                
                [f, x] = ksdensity(yLengths);
                [~, maxIndex] = max(f);
                medianGlintVector(2) = x(maxIndex);
                
                
            end
        else
            unsortedGlintVector = [];
            medianGlintVector = p.Results.medianGlintVector;
        end
        % sort tracked glints according to the glint vector main direction:
        % make main direction component positive by swapping the tracked
        % centroid order if necessary
        switch p.Results.glintsMainDirection
            case {'y', 'both'}

                for ii = 1:size(unsortedGlintVector,1)

                    
                    [lengthError, candidateGlintVectorIdx] = min(sqrt((abs(unsortedGlintVector(ii,2)) - medianGlintVector(2)).^2 + (abs(unsortedGlintVector(ii,1)) - medianGlintVector(1)).^2));
                    
                    desiredLength = sqrt(medianGlintVector(2)^2 + medianGlintVector(1)^2);
                     
                    if isempty(p.Results.threshold) || (~isempty(p.Results.threshold) && p.Results.threshold > lengthError/desiredLength)
                        
                        if unsortedGlintVector(ii,2)>0
                            glintData_X(framesWithExpectedCentroids(ii),1) = centroidsByFrame_X(framesWithExpectedCentroids(ii),1);
                            glintData_X(framesWithExpectedCentroids(ii),2) = centroidsByFrame_X(framesWithExpectedCentroids(ii),2);
                            
                            glintData_Y(framesWithExpectedCentroids(ii),1) = centroidsByFrame_Y(framesWithExpectedCentroids(ii),1);
                            glintData_Y(framesWithExpectedCentroids(ii),2) = centroidsByFrame_Y(framesWithExpectedCentroids(ii),2);
                            
                        else
                            glintData_X(framesWithExpectedCentroids(ii),1) = centroidsByFrame_X(framesWithExpectedCentroids(ii),2);
                            glintData_X(framesWithExpectedCentroids(ii),2) = centroidsByFrame_X(framesWithExpectedCentroids(ii),1);
                            
                            glintData_Y(framesWithExpectedCentroids(ii),1) = centroidsByFrame_Y(framesWithExpectedCentroids(ii),2);
                            glintData_Y(framesWithExpectedCentroids(ii),2) = centroidsByFrame_Y(framesWithExpectedCentroids(ii),1);
                        end
                    end
                end
                
            case 'x'
                for ii = 1:length(unsortedGlintVector)
                    if unsortedGlintVector(ii,1)>0
                        glintData_X(framesWithExpectedCentroids(ii),1) = centroidsByFrame_X(framesWithExpectedCentroids(ii),1);
                        glintData_X(framesWithExpectedCentroids(ii),2) = centroidsByFrame_X(framesWithExpectedCentroids(ii),2);
                        
                        glintData_Y(framesWithExpectedCentroids(ii),1) = centroidsByFrame_Y(framesWithExpectedCentroids(ii),1);
                        glintData_Y(framesWithExpectedCentroids(ii),2) = centroidsByFrame_Y(framesWithExpectedCentroids(ii),2);
                        
                    else
                        glintData_X(framesWithExpectedCentroids(ii),1) = centroidsByFrame_X(framesWithExpectedCentroids(ii),2);
                        glintData_X(framesWithExpectedCentroids(ii),2) = centroidsByFrame_X(framesWithExpectedCentroids(ii),1);
                        
                        glintData_Y(framesWithExpectedCentroids(ii),1) = centroidsByFrame_Y(framesWithExpectedCentroids(ii),2);
                        glintData_Y(framesWithExpectedCentroids(ii),2) = centroidsByFrame_Y(framesWithExpectedCentroids(ii),1);
                    end
                end
        end
        
    otherwise
        error('Sorry, unable to track this many glints at the moment.')
end

% now, find frames with more centroids than expected
framesWithMoreCentroids = find (centroidsInEachFrame>p.Results.numberOfGlints);

% select the centroids that are most likely to be glints
if ~isempty(framesWithMoreCentroids)
    
    switch p.Results.numberOfGlints
        case 1
            % If we have a perimeter, pick the glint that is closest to the
            % center of the pupil. If not, we use the median position of
            % the centroids.
            if ~isempty(p.Results.perimeterFileName)
                load(p.Results.perimeterFileName,'perimeter');
                target_X = cellfun(@(x) mean(x.Xp),perimeter.data);
                target_Y = cellfun(@(x) mean(x.Yp),perimeter.data);
            else
                target_X = repmat(nanmedian(centroidsByFrame_X(:)),1,nFrames);
                target_Y = repmat(nanmedian(centroidsByFrame_Y(:)),1,nFrames);
            end
            
            % loop through the frames with too many centroids and check
            % which value is closer to the target.
            for ii = 1:length(framesWithMoreCentroids)
                % find the centroid closest to the target
                [~,glintIDX] = min(sqrt(...
                    (centroidsByFrame_X(framesWithMoreCentroids(ii),:) - target_X(framesWithMoreCentroids(ii))).^2 + ...
                    (centroidsByFrame_Y(framesWithMoreCentroids(ii),:)- target_Y(framesWithMoreCentroids(ii))).^2));
                % store values for that centroid as glint
                if glintIDX>=1 && glintIDX<=p.Results.centroidsAllocation
                    glintData_X(framesWithMoreCentroids(ii)) = centroidsByFrame_X(framesWithMoreCentroids(ii),glintIDX);
                    glintData_Y(framesWithMoreCentroids(ii)) = centroidsByFrame_Y(framesWithMoreCentroids(ii),glintIDX);
                end
            end
            
        case 2
            % loop through frames
            for ii = 1:length(framesWithMoreCentroids)                
                % get all possible combinations of glint vector main
                % length, and select the couple of glints that generate a
                % vector main length closest to the average.
                switch p.Results.glintsMainDirection
                    case 'y'
                        theseCentroids = (centroidsByFrame_Y(framesWithMoreCentroids(ii),:));
                        
                    case 'x'
                        theseCentroids = (centroidsByFrame_X(framesWithMoreCentroids(ii),:));
                        
                    case 'both'
                        theseCentroids(1,:) = (centroidsByFrame_Y(framesWithMoreCentroids(ii),:));
                        theseCentroids(2,:) = (centroidsByFrame_X(framesWithMoreCentroids(ii),:));
                    otherwise
                        error ('Main direction must be ''x'' or ''y'' for 2 glints case')
                end
                switch p.Results.glintsMainDirection
                    case {'x' 'y'}
                        % compute all possible main lengths
                        mainLengths = pdist(theseCentroids',@(x,y) x-y);
                        
                        % find main length closest to the average
                        [lengthError, candidateGlintVectorIdx] = min(abs(abs(mainLengths) - medianGlintVector(2)));
                        desiredLength = medianGlintVector(2);
                        
                        
                    case 'both'

                        mainLengthsY = pdist(theseCentroids(1,:)');
                        mainLengthsX = pdist(theseCentroids(2,:)');
                        
                        [lengthError, candidateGlintVectorIdx] = min(sqrt((abs(mainLengthsY) - medianGlintVector(2)).^2 + (abs(mainLengthsX) - medianGlintVector(1)).^2));
                        
                        desiredLength = sqrt(medianGlintVector(2)^2 + medianGlintVector(1)^2);
                        % find main length closest to the average
                        %[minimumLength, candidateGlintVectorIdx] = min(abs(abs(mainLengths) - (sqrt(medianGlintVector(2)^2 + medianGlintVector(1)^2))));
                        
                end
                
                
                if isempty(p.Results.threshold) || (~isempty(p.Results.threshold) && p.Results.threshold > lengthError/desiredLength)
                    % get the indexes of the appropriate centroids
                    if candidateGlintVectorIdx < length(theseCentroids)
                        glintsTemp_Y = [(centroidsByFrame_Y(framesWithMoreCentroids(ii),1)) (centroidsByFrame_Y(framesWithMoreCentroids(ii),(candidateGlintVectorIdx+1)))];
                        glintsTemp_X = [(centroidsByFrame_X(framesWithMoreCentroids(ii),1)) (centroidsByFrame_X(framesWithMoreCentroids(ii),(candidateGlintVectorIdx+1)))];
                        
                    elseif p.Results.centroidsAllocation == 5 && (ismember(candidateGlintVectorIdx,[5 6 7]))
                        glintsTemp_Y = [(centroidsByFrame_Y(framesWithMoreCentroids(ii),2)) (centroidsByFrame_Y(framesWithMoreCentroids(ii),(candidateGlintVectorIdx -2)))];
                        glintsTemp_X = [(centroidsByFrame_X(framesWithMoreCentroids(ii),2)) (centroidsByFrame_X(framesWithMoreCentroids(ii),(candidateGlintVectorIdx -2)))];
                        
                    elseif p.Results.centroidsAllocation == 5 && (ismember(candidateGlintVectorIdx,[8 9]))
                        glintsTemp_Y = [(centroidsByFrame_Y(framesWithMoreCentroids(ii),3)) (centroidsByFrame_Y(framesWithMoreCentroids(ii),(candidateGlintVectorIdx -4)))];
                        glintsTemp_X = [(centroidsByFrame_X(framesWithMoreCentroids(ii),3)) (centroidsByFrame_X(framesWithMoreCentroids(ii),(candidateGlintVectorIdx -4)))];
                        
                    elseif p.Results.centroidsAllocation == 5 && candidateGlintVectorIdx == 10
                        glintsTemp_Y = [(centroidsByFrame_Y(framesWithMoreCentroids(ii),4)) (centroidsByFrame_Y(framesWithMoreCentroids(ii),5))];
                        glintsTemp_X = [(centroidsByFrame_X(framesWithMoreCentroids(ii),4)) (centroidsByFrame_X(framesWithMoreCentroids(ii),5))];
                        
                    else
                        error ('Option currently unavailable, please set centroidsAllocation = 5 and run again')
                    end
                    
                    % sort the glints
                    
                    switch p.Results.glintsMainDirection
                        case {'y', 'both'}
                            glintsTemp = glintsTemp_Y;
                            
                        case 'x'
                            glintsTemp = glintsTemp_X;
                            
                        otherwise
                            error ('Main direction must be ''x'' or ''y'' for 2 glints case')
                    end
                    
                    if diff(glintsTemp) >0
                        glintData_X(framesWithMoreCentroids(ii),1) = glintsTemp_X(2);
                        glintData_X(framesWithMoreCentroids(ii),2) = glintsTemp_X(1);
                        
                        glintData_Y(framesWithMoreCentroids(ii),1) = glintsTemp_Y(2);
                        glintData_Y(framesWithMoreCentroids(ii),2) = glintsTemp_Y(1);
                    else
                        glintData_X(framesWithMoreCentroids(ii),2) = glintsTemp_X(2);
                        glintData_X(framesWithMoreCentroids(ii),1) = glintsTemp_X(1);
                        
                        glintData_Y(framesWithMoreCentroids(ii),2) = glintsTemp_Y(2);
                        glintData_Y(framesWithMoreCentroids(ii),1) = glintsTemp_Y(1);
                    end
                    
                    clear glintsTemp glintsTemp_Y glintsTemp_X theseCentroids
                end
            end
            
        otherwise
            error('Sorry, unable to track this many glints at the moment.')
    end %switch number of glints
end

% final glint clean up. If we found a glints in one frame, but not in
% either around it, assume those isolated glints are spurious and should be
% removed

if p.Results.removeIsolatedGlints
    for ii = 2:length(glintData_X) - 1
        if sum(isnan(glintData_X(ii-1,:))) == 2 && sum(isnan(glintData_X(ii+1,:))) == 2
            glintData_X(ii,1) = nan;
            glintData_X(ii,2) = nan;

            glintData_Y(ii,1) = nan;
            glintData_Y(ii,2) = nan;

        end
        
    end
end

% finally, in the case there are fewer centroids than expected, we assume
% that the real glints are not visible in the frame (e.g. it is a blink
% frame) and the centroids are picked up randomly from bright spots on the
% image. We therefore leave those frames with NaN values for the glint
% locations.


%% save out all data in glintData struct
glintData.X = glintData_X;
glintData.Y = glintData_Y;
glintData.meta = p.Results;
glintData.meta.centroidsByFrame.X = centroidsByFrame_X;
glintData.meta.centroidsByFrame.Y = centroidsByFrame_Y;
glintData.meta.coordinateSystem = 'intrinsicCoordinates(pixels)';


% close the video object
clear videoInObj

% save out a mat file with the glint tracking data
if ~p.Results.displayMode
    save (glintFileName, 'glintData')
else
    close(figureHandle);
end

% report completion of analysis
if p.Results.verbose
    toc
    fprintf('\n');
end


end % main function

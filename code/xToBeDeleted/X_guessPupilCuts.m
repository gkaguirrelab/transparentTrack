function framesToCut = guessPupilCuts(perimeterVideoName,glintFileName,blinkFrames,varargin)

% this function makes the best guess on which frames need cutting of part
% of the pupil perimeter for a better ellipse fit. The guess is based on a
% simple attempt of fitting an ellipse to the perimeter. If the distance
% between the points of the perimeter and those of the ellipse (i.e. the
% fitting error) exceeeds an arbitrary threshold, the routine will try
% cutting out portions of the perimeter and return the cut that yelds the
% following information: which frames to cut, what was their original
% fitting error, which is the suggested cut, what is the fitting error with
% the suggested cut.
%
% Note that the cut depends on the glint location.
%
%
% Output
% ======
%       framesToCut = [frameNumber originalError Ucut Rcut newError]
%
% Input params
% ============
%       perimeterVideoName : path to the perimeter video file to track
%       glintFileName : path to the matFile of glint location.
%       blinkFrames : index of blink frames
%
% Options
% =======
%       errorThreshold : the distance error tolerated before attempting to
%           cut (default: 10)
%       TO BE DEVELOPED availableCuts : which set of cuts to include when trying to cut the
%           perimeter. (default: 'all')
%
%
%
% Usage example
% =============
%  framesToCut = guessPupilCuts(perimeterVideoName,glintFileName,blinkFrames);



%% parse input and define variables

p = inputParser;
% required input
p.addRequired('perimeterVideoName',@isstr);
p.addRequired('glintFileName',@isstr);
p.addRequired('blinkFrames',@isnumeric);

% optional inputs
errorThresholdDefault = 10;
p.addParameter('errorThreshold', errorThresholdDefault, @isnumeric);

%parse
p.parse(perimeterVideoName, glintFileName, blinkFrames, varargin{:})

% define optional variables values
errorThreshold = p.Results.errorThreshold;

%% load and prepare data

% load glint
load(glintFileName)

% pupilPerimeter
inObj = VideoReader(perimeterVideoName);

% get number of frames
numFrames = floor(inObj.Duration*inObj.FrameRate);

% NOTE function assumes that blinkFrames is a numeric variable for now.

%% initialize framesToCut

% to begin with, we initialize frameToCut to have as many rows as the
% number of frames. Before saving out the output, we will trim the unused row
% out.

framesToCut = nan(numFrames, 5);


%% loop through all frames the first time. Try to apply a horizontal cut for high error frames.

% initialize progress bar
progBar = ProgressBar(numFrames,'Guessing pupil cuts 1/2...');

for ii = 1:numFrames
    
    % skip the frames marked as blinks
    if ismember(ii,blinkFrames)
        if ~mod(ii,10);progBar(ii);end % update progressbar
        continue
    end
    
    % read the frame
    thisFrame = read(inObj,ii);
    % NOTE: in fugure releases this will need to be replaced. However,
    % readFrame currently only reads frames progressively, so we keep using
    % read instead to skip the frames we don't care about.
    
    % convert frame to binary image and index the perimeter points
    thisFrame = rgb2gray (thisFrame);
    binP = imbinarize(thisFrame);
    [Yp, Xp] = ind2sub(size(binP),find(binP));
    
    % check that the frame is not empty
    if isempty (Yp)
         if ~mod(ii,10);progBar(ii);end % update progressbar
        continue
    end
    
    % fit an ellipse to the full perimeter using the quadFit toolbox
    try
        Epi = ellipsefit_direct(Xp,Yp);
        [~,d,~,~] = ellipse_distance(Xp, Yp, Epi);
        originalFittingError = nanmedian(sqrt(sum(d.^2)));
    catch ME
    end
    if  exist ('ME', 'var')
        framesToCut(ii,:) = [ii 0 0 0 0];  % the all zero output means: fitting error.
        clear ME
        if ~mod(ii,10);progBar(ii);end % update progressbar
        continue
    end
    
    % check the fitting error
    if originalFittingError <= errorThreshold
        if ~mod(ii,10);progBar(ii);end % update progressbar
        continue
    else
        % try cuts in a prioritized order until you find the best one
        % (prioritize cutting the least amount of pixels).
        
        % get glint position for this frame
        Xg = glintData.X(ii);
        Yg = glintData.Y(ii);
        
        % FIRST WE TRY HORIZONTAL CUTS (most common occurrence)
        % define vertical cuts
        U = flip(0:1:round(Yg - (min(Yp)+1)));
        R = Inf;
        
        % initialize values for while loop
        cc = 1;
        newFittingError = errorThreshold;
        
        % while loop to find the best cut
        while newFittingError >= errorThreshold && cc <= length(U)
            [binPcut] = cutPupil (binP,U(cc),R,Xg,Yg);
            [Yc, Xc] = ind2sub(size(binPcut),find(binPcut));
            
            try
                Epi = ellipsefit_direct(Xc,Yc);
                [~,d,~,~] = ellipse_distance(Xc, Yc, Epi);
                newFittingError = nanmedian(sqrt(sum(d.^2)));
            catch ME
            end
            if  exist ('ME', 'var')
                framesToCut(ii,:) = [ii originalFittingError 0 0 0];  % the all zero output means: fitting error.
                clear ME
                continue
            end
            allErrors(cc) = newFittingError;
            cc = cc+1;
        end
        if cc == length(U) + 1
            % take best cut
            try
                [bestError, beIdx] = min(allErrors);
                framesToCut(ii,:) = [ii originalFittingError U(beIdx) R bestError];
            catch ME
            end
            if  exist ('ME', 'var')
                framesToCut(ii,:) = [ii originalFittingError 0 0 0];  % the all zero output means: fitting error.
                clear ME
                continue
            end
        else
            framesToCut(ii,:) = [ii originalFittingError U(cc) R newFittingError];
        end
    end
    if ~mod(ii,10);progBar(ii);end % update progressbar
end % loop through all frames

% remove NaN rows
framesToCut(any(isnan(framesToCut),2),:) = [];

%% NOW WE TRY HORIZONTAL CUTS

% we just look of what is still bad after the vertical cuts
fittingErrorIdx = find(framesToCut(:,5) == 0);
highErrorIdx = find(framesToCut(:,5)> errorThreshold);

highErrorIdx = union(fittingErrorIdx,highErrorIdx);
if ~isempty(highErrorIdx)
    
    highErrorFrames = framesToCut(highErrorIdx,1);
    
    % loop in highErrorFrames and do vertical cuts
    
    progBar = ProgressBar(highErrorFrames,'Guessing pupil cuts 2/2...');
    
    for ee = 1: length(highErrorFrames)
        
        % read the frame
        thisFrame = read(inObj,highErrorFrames(ee));
        % convert frame to binary image and index the perimeter points
        thisFrame = rgb2gray (thisFrame);
        binP = imbinarize(thisFrame);
        [Yp, Xp] = ind2sub(size(binP),find(binP));
        
        % get glint position for this frame
        Xg = glintData.X(highErrorFrames(ee));
        Yg = glintData.Y(highErrorFrames(ee));
        
        % define horizontal cuts
        R = flip(0:1:round((max(Xp)-1) - Xg ));
        U = Inf;
        
        % initialize for while loop
        cc = 1;
        newFittingError = errorThreshold;
        
        while newFittingError >= errorThreshold && cc <= length(R)
            [binPcut] = cutPupil (binP,U,R(cc),Xg,Yg);
            [Yc, Xc] = ind2sub(size(binPcut),find(binPcut));
            
            try
                Epi = ellipsefit_direct(Xc,Yc);
                [~,d,~,~] = ellipse_distance(Xc, Yc, Epi);
                newFittingError = nanmedian(sqrt(sum(d.^2)));
            catch ME
            end
            if  exist ('ME', 'var')
                newFramesToCut(ee,:) = [highErrorFrames(ee) 0 0 0];  % the all zero output means: fitting error.
                clear ME
                 if ~mod(ee,10);progBar(ee);end % update progressbar
                continue
            end
            allErrors(cc) = newFittingError;
            cc = cc+1;
        end
        if cc == length(R) + 1
            % take best cut
            try
                [bestError, beIdx] = min(allErrors);
                newFramesToCut(ee,:) = [highErrorFrames(ee) U R(beIdx) bestError];
            catch ME
            end
            if  exist ('ME', 'var')
                newFramesToCut(ee,:) = [highErrorFrames(ee) 0 0 0];  % the all zero output means: fitting error.
                clear ME
                 if ~mod(ee,10);progBar(ee);end % update progressbar
                continue
            end
        else
            newFramesToCut(ee,:) = [highErrorFrames(ee) U R(cc) newFittingError];
        end
       if ~mod(ee,10);progBar(ee);end % update progressbar  
    end
    
    % now compare the error from the horizontal cut to this one. If
    % this one is better,take this cut instead.
    
    for kk = 1:size(newFramesToCut,1)
        if framesToCut(highErrorIdx(kk),5) > newFramesToCut(kk,4) && newFramesToCut(kk,4) ~= 0
            framesToCut(highErrorIdx(kk),3) = newFramesToCut(kk,2);
            framesToCut(highErrorIdx(kk),4) = newFramesToCut(kk,3);
            framesToCut(highErrorIdx(kk),5) = newFramesToCut(kk,4);
        end
    end     
end

%% At this point: 
% all that's left is to find the frames that still have a high error or had
% a fitting error (all zero in the control file). We can either apply a
% diagonal cut with an iterative procedure, or have a human select
% appropriate diagonal cut for them.


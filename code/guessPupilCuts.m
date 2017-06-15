function framesToCut = guessPupilCuts(perimeterVideo,glintFile,blinkFrames,varargin)

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
%       perimeterVideo : path to the perimeter video file to track
%       glintFile : path to the matFile of glint location.
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
%  framesToCut = guessPupilCuts(perimeterVideo,glintFile,blinkFrames);



%% parse input and define variables

p = inputParser;
% required input
p.addRequired('perimeterVideo',@isstr);
p.addRequired('glintFile',@isstr);
p.addRequired('blinkFrames',@isnumeric);

% optional inputs
errorThresholdDefault = 10;
p.addParameter('errorThreshold', errorThresholdDefault, @isnumeric);

%parse
p.parse(perimeterVideo, glintFile, blinkFrames, varargin{:})

% define optional variables values
errorThreshold = p.Results.errorThreshold;

%% load and prepare data

% load glint
load(glintFile)

% pupilPerimeter
inObj = VideoReader(perimeterVideo);

% get number of frames
numFrames = floor(inObj.Duration*inObj.FrameRate);

% NOTE function assumes that blinkFrames is a numeric variable for now.

%% initialize framesToCut

% to begin with, we initialize frameToCut to have as many rows as the
% number of frames. Before saving out the output, we will trim the unused row
% out.

framesToCut = nan(numFrames, 5);

%% guess the cuts

for ii = 1:numFrames
    
    % skip the frames marked as blinks
    if ismember(ii,blinkFrames)
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
        continue
    end
    
    % check the fitting error
    if originalFittingError <= errorThreshold
        continue
    else
        % try cuts in a prioritized order until you find the best one
        % (prioritize cutting the least amount of pixels).
        
        % get glint position for this frame
        Xg = glint.X(ii);
        Yg = glint.Y(ii);
        
        % define vertical cuts
        U = flip(0:1:round(Yg - (min(Yp)+1)));
        R = Inf;
        
        % initialize for while loop
        cc = 1;
        newFittingError = errorThreshold;
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
    
end

% remove NaN rows
framesToCut(any(isnan(framesToCut),2),:) = [];







function [frameSet, gazeTargets] = gazePost(videoStemName, varargin)
% Identify a fixation frame that can be used to sync sceneGeometry
%
% Syntax:
%  [frameSet, gazeTargets] = selectFrames.gazePost(videoStemName)
%
% Description:
%   Positioning an eye model in a scene requires the selection of
%   informative frames of the acquisition to guide the alignment.
%
%   % For some acqusitions (i.e., retinotopy), the subject was asked to
%   stare at a fixation point in the center of the screen after the start
%   of the acquisition. This could also work for movie viewing, if we are
%   willing to assume that the median gaze position during the first (e.g.)
%   10 seconds of a movie is the center of the screen. This routine finds a
%   frame from the pre-scan acquisition that likely has a gazeTarget values
%   of [0; 0].
%
% Inputs:
%	videoStemName         - Char vector. Full path to video file from which
%                           the scene observations have been derived. The
%                           stem name should omit the "_gray.avi" suffix
%                           that is usually present in the names of these
%                           video files.
%
% Optional key-value pairs:
%  'eyePositionTargetLengthFrames' - Scalar. The number of sequential 
%                           frames from the target acquisition that will be
%                           found and used to define the position of the
%                           eye to be fit.
%  'gazeErrorThreshTol'   - Scalar. The run of frames must have a deviation
%                           of less than this value. The precise meaning of
%                           the value will differ for the different
%                           alignment methods.
%
% Outputs:
%   frameSet              - Scalar that specifies a frame index
%                           (indexed from 1).
%   gazeTargets           - A 2x1 matrix that provides the positions, in
%                           degrees of visual angle, of the likely fixation
%                           position [0;0] of the eye for this frame.
%


%% input parser
p = inputParser; p.KeepUnmatched = true;

% Required
p.addRequired('videoStemName',@ischar);

p.addParameter('eyePositionTargetLengthFrames',30,@isscalar);
p.addParameter('gazeErrorThreshTol',0.25,@isscalar);

% parse
p.parse(videoStemName, varargin{:})

% Load the timebase, pupilData, perimeter, and relative camera position for
% this acquisition
load([videoStemName '_timebase.mat'],'timebase');
load([videoStemName '_pupil.mat'],'pupilData');
load([videoStemName '_glint.mat'],'glintData');

% Identify the acqStartTimeFixed, which is the time point at which the
% fMRI acquisition began
[~, acqStartFrameFixed] = min(abs(timebase.values));


% For some acqusitions (i.e., retinotopy), the subject was asked to
% stare at a fixation point in the center of the screen after the
% start of the acquisition. This could also work for movie viewing,
% if we are willing to assume that the median gaze position during
% the first (e.g.) 10 seconds of a movie is the center of the screen.
% Find the period after to the start of the scan when the eye was
% in the most consistent position, and closest to the median
% position.
windowStart = acqStartFrameFixed;
windowEnd = acqStartFrameFixed+600;
gazeX = pupilData.initial.ellipses.values(windowStart:windowEnd,1);
gazeY = pupilData.initial.ellipses.values(windowStart:windowEnd,2);
medianX = nanmedian(gazeX);
medianY = nanmedian(gazeY);

% Assuming that the eye had a central tendency of fixation upon the
% center of the screen, this vector expresses the deviation of
% fixation on any given frame from the screen center.
eyeMatchError = sqrt(sum([gazeX-medianX; gazeY-medianY].^2,2));

% Set any period without a glint to have an arbitrarily high error
eyeMatchError(isnan(glintData.X)) = 1e20;

% Find the minimum fixation error threshold that results in a run of
% consecutive frames at fixation of the target length.
targetLength = p.Results.eyePositionTargetLengthFrames;

% Anonymous function to grab the indicies of when runs of frames begin
pullStartIndices = @(vec) vec(1:2:end-1);

% Anonymous function to grab the length of each run
pullRunLengths = @(vec) vec(2:2:end)-pullStartIndices(vec);

% Anonynous function that provides the lengths of runs of frames for which
% the eyeMatch error is below a threshold.
runStarts = @(thresh) find(diff([0,(eyeMatchError < thresh)',0]==1));

% Set the fzero search options
options = optimset('fzero');
options.Display = 'off';

% The starting point for the search
x0 = 0.5;

% Adjust the targetLength as needed to achieve a threshVal below threshTol.
stillWorking = true;
while stillWorking
    % An objective function that expresses the difference of the longest
    % run length from the target run length (e.g., 30 frames long). The
    % business with the min([1e6 ...]) is to handle the case when the run
    % set is empty, and thus would otherwise return an empty variable for
    % the objective.
    myObj = @(thresh) min([1e6, (targetLength - max(pullRunLengths(runStarts(thresh))))]);
    
    % Perform the search
    threshValFixed = fzero(myObj,x0,options);
    
    % Check to see if the search has met our threshold criteria, or if we
    % have run the targetLength down as short as it can go.
    if threshValFixed < p.Results.gazeErrorThreshTol || targetLength == 1
        stillWorking = false;
    else
        targetLength = targetLength-1;
    end
end

% Check if we found a solution
if ~isfinite(threshValFixed)
    error('gazePost:noGoodFrame','Unable to find a suitable set of frames from the acquisition')
end

% Find the start point of this run of frames
runLengths = pullRunLengths(runStarts(threshValFixed));
runIndices = pullStartIndices(runStarts(threshValFixed))+windowStart-1;
runLengthFixed = targetLength-myObj(threshValFixed);
startIndexFixed = runIndices(runLengths == runLengthFixed);
startIndexFixed = startIndexFixed(1);
frameSetFixed = startIndexFixed:startIndexFixed+runLengthFixed-1;

% Grab RMSE of the ellipse fit to each frame
ellipseRMSE = pupilData.initial.ellipses.RMSE;

% Set any frame without a glint to nan
ellipseRMSE(isnan(glintData.X)) = nan;

% Work with just the frameSet
ellipseRMSE = ellipseRMSE(frameSetFixed);

% Check if we found a solution
if all(isnan(ellipseRMSE))
    error('gazePre:noGlint','There is no glint available during the gazePre period')
end

% Find the frame with the lowest ellipse RMSE during the target window
frameSet = startIndexFixed + find(ellipseRMSE == min(ellipseRMSE)) - 1;

% The gaze target is presumed to be the center of the screen [0 0]
gazeTargets = [0; 0];



end



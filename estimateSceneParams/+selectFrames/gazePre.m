function [frameSet, gazeTargets, relativeCameraPosition] = gazePre(videoStemName, varargin)



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
load([videoStemName '_relativeCameraPosition.mat'],'relativeCameraPosition');

% Identify the acqStartTimeFixed, which is the time point at which the
% fMRI acquisition began
[~, acqStartFrameFixed] = min(abs(timebase.values));

        % For some acquisitions, the subject was asked to stare at a
        % fixation point in the center of the screen prior to the start of
        % the acquisition. Find the period prior to the start of the scan
        % when the eye was in the most consistent position, and closest to
        % the median position.
        windowStart = 1;
        windowEnd = acqStartFrameFixed;
        gazeX = pupilData.initial.ellipses.values(windowStart:windowEnd,1);
        gazeY = pupilData.initial.ellipses.values(windowStart:windowEnd,2);
        medianX = nanmedian(gazeX);
        medianY = nanmedian(gazeY);
        
        % Assuming that the eye had a central tendency of fixation upon the
        % center of the screen, this vector expresses the deviation of
        % fixation on any given frame from the screen center.
        eyeMatchError = sqrt(sum([gazeX-medianX; gazeY-medianY].^2,2));
        x0 = 0.5;
        
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
    warning('Unable to find a suitable set of frames from the acquisition')
    return
end

% Find the start point of this run of frames
runLengths = pullRunLengths(runStarts(threshValFixed));
runIndices = pullStartIndices(runStarts(threshValFixed))+windowStart-1;
runLengthFixed = targetLength-myObj(threshValFixed);
startIndexFixed = runIndices(runLengths == runLengthFixed);
startIndexFixed = startIndexFixed(1);
frameSetFixed = startIndexFixed:startIndexFixed+runLengthFixed-1;


% Find the frame with the lowest ellipse RMSE during the target window
ellipseRMSE = pupilData.initial.ellipses.RMSE(frameSetFixed);
frameSet = startIndexFixed + find(ellipseRMSE == min(ellipseRMSE)) - 1;

% The gaze target is presumed to be the center of the screen [0 0]
gazeTargets = [0; 0];

relativeCameraPosition = relativeCameraPosition.values(:,frameSet);

end
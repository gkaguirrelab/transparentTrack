function prepareGazeCalibrationData (LTdatFileName,gazeDataFileName,varargin)

% prepareGazeCalibrationData

% HEADER


%% input parser









%% load data from ltDatFile
% get target location
calParams.targets.X     = LTdata.targets(:,1); % mm on screen, screen center = 0
calParams.targets.Y     = LTdata.targets(:,2); % mm on screen, screen center = 0

% if the livetrack failed to track one target, it will be a NaN. Replace
% the NaN with the missing target location. NOTE THAT THIS ASSUMES THAT
% ONLY ONE TARGET IS MISSING

% find nan target
nanIDX = find(isnan(calParams.targets.X));

if ~isempty(nanIDX)
   if length(nanIDX) == 2
       warning ('There is more than one NaN target! Calibration might fail.')
   end
    % available targets locations
    high = max(LTdata.targets(:,1));
    center = 0;
    low = min(LTdata.targets(:,1));
    
    allLocations = [...
        high high; ...
        high center; ...
        high low; ...
        center high; ...
        center center; ...
        center low; ...
        low high; ...
        low center; ...
        low low; ...
        ];
    
    % find value of the nan target (is the one missing when comparing the
    % targets to allLocations)

%         allLocations = flipud(allLocations);
    
    missingTarget = find(~ismember(allLocations,[calParams.targets.X calParams.targets.Y], 'rows'));
    
    % replace NaN target value with missing value
    
    calParams.targets.X(nanIDX) = allLocations(missingTarget, 1);
    calParams.targets.Y(nanIDX) = allLocations(missingTarget, 2);
end

% get viewing distance


% get dot times (if it was saved)



%% Switch cases
% got raw data and dot times

% FROM OLD FUNCTION
% get each target duration on screen
TarTimesFromStart  = (round(LTdata.dotTimes - rawVidStart) * fps - 40);

targetDurSec = diff(LTdata.dotTimes); % target duration in seconds

targetDurFrames = round(targetDurSec * fps);
videoStartFrames = round(LTdata.dotTimes(1) - rawVidStart) * fps;
totalFrames = round(LTdata.dotTimes(end) - rawVidStart) * fps;


% % remove some of the boundaries
removeFrames = round(targetDurFrames ./10);

% make moving mean on window with 1/3 of samples
window = round(targetDurFrames ./3);
for ct = 1 :9
    movMeanP.X= movmedian(trackData.pupil.X(TarTimesFromStart(ct) + removeFrames(ct) :TarTimesFromStart(ct+1) - removeFrames(ct)), window(ct),'omitnan','Endpoints','discard');
    movMeanP.Y = movmedian(trackData.pupil.Y(TarTimesFromStart(ct) + removeFrames(ct) :TarTimesFromStart(ct+1) - removeFrames(ct)), window(ct),'omitnan','Endpoints','discard');
    movMeanG.X = movmedian(trackData.glint.X(TarTimesFromStart(ct) + removeFrames(ct) :TarTimesFromStart(ct+1) - removeFrames(ct)), window(ct),'omitnan','Endpoints','discard');
    movMeanG.Y = movmedian(trackData.glint.Y(TarTimesFromStart(ct) + removeFrames(ct) :TarTimesFromStart(ct+1) - removeFrames(ct)), window(ct),'omitnan','Endpoints','discard');
    
    movStdP.X = movstd(trackData.pupil.X(TarTimesFromStart(ct) + removeFrames(ct) :TarTimesFromStart(ct+1) - removeFrames(ct)), window(ct),'omitnan','Endpoints','discard');
    movStdP.Y = movstd(trackData.pupil.Y(TarTimesFromStart(ct) + removeFrames(ct) :TarTimesFromStart(ct+1) - removeFrames(ct)), window(ct),'omitnan','Endpoints','discard');
    movStdG.X = movstd(trackData.glint.X(TarTimesFromStart(ct) + removeFrames(ct) :TarTimesFromStart(ct+1) - removeFrames(ct)), window(ct),'omitnan','Endpoints','discard');
    movStdG.Y = movstd(trackData.glint.Y(TarTimesFromStart(ct) + removeFrames(ct) :TarTimesFromStart(ct+1) - removeFrames(ct)), window(ct),'omitnan','Endpoints','discard');
    
    % ger minimum std
    [minPX,Xidx] = min(movStdP.X);
    [minPY,Yidx] = min(movStdP.Y);
    
    % take according to x
    meanP.X(ct) = movMeanP.X(Xidx);
    meanP.Y(ct) = movMeanP.Y(Xidx);
    meanG.X(ct) = movMeanG.X(Xidx);
    meanG.Y(ct) = movMeanG.Y(Xidx);
end




% got raw data no dot times




% no raw data at all
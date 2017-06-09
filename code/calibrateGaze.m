
function calParams = calibrateGaze (dropboxDir, params)


if ~isfield(params,'viewDist')
    params.viewDist = 1065;
end


fps = 60; % fps


trackData = load(params.trackedData);% = load ('/Volumes/Bay_2_data/giulia/Dropbox-Aguirre-Brainard-Lab/TOME_processing/session2_spatialStimuli/TOME_3008/103116/EyeTracking/GazeCal02_testTracking.mat');
LTdata = load(params.LTcal); %load('/Volumes/Bay_2_data/giulia/Dropbox-Aguirre-Brainard-Lab/TOME_data/session2_spatialStimuli/TOME_3008/103116/EyeTracking/GazeCal02_LTdat.mat');

load(params.rawVidStart) %'/Volumes/Bay_2_data/giulia/Dropbox-Aguirre-Brainard-Lab/TOME_data/session2_spatialStimuli/TOME_3008/103116/EyeTracking/GazeCal02_rawVidStart.mat');

%% get the fixation data

% get each target duration on screen
TarTimesFromStart  = (round(LTdata.dotTimes - rawVidStart) * fps - 40);

targetDurSec = diff(LTdata.dotTimes); % target duration in seconds

targetDurFrames = round(targetDurSec * fps);
videoStartFrames = round(LTdata.dotTimes(1) - rawVidStart) * fps;
totalFrames = round(LTdata.dotTimes(end) - rawVidStart) * fps;


% % remove some of the boundaries
removeFrames = round(targetDurFrames ./10);
%
% for ct = 1 :9
%     meanP.X(ct) = nanmedian(trackData.pupil.X(TarTimesFromStart(ct) + removeFrames(ct) :TarTimesFromStart(ct+1) - removeFrames(ct)));
%     meanP.Y(ct) = nanmedian(trackData.pupil.Y(TarTimesFromStart(ct) + removeFrames(ct) :TarTimesFromStart(ct+1) - removeFrames(ct)));
%     meanG.X(ct) = nanmedian(trackData.glint.X(TarTimesFromStart(ct) + removeFrames(ct) :TarTimesFromStart(ct+1) - removeFrames(ct)));
%     meanG.Y(ct) = nanmedian(trackData.glint.Y(TarTimesFromStart(ct) + removeFrames(ct) :TarTimesFromStart(ct+1) - removeFrames(ct)));
%
%     stdP.X(ct) = nanstd(trackData.pupil.X(TarTimesFromStart(ct) + removeFrames(ct) :TarTimesFromStart(ct+1) - removeFrames(ct)));
%     stdP.Y(ct) = nanstd(trackData.pupil.Y(TarTimesFromStart(ct) + removeFrames(ct) :TarTimesFromStart(ct+1) - removeFrames(ct)));
%     stdG.X(ct) = nanstd(trackData.glint.X(TarTimesFromStart(ct) + removeFrames(ct) :TarTimesFromStart(ct+1) - removeFrames(ct)));
%     stdG.Y(ct) = nanstd(trackData.glint.Y(TarTimesFromStart(ct) + removeFrames(ct) :TarTimesFromStart(ct+1) - removeFrames(ct)));
%
% end

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


% get data ready for calibration
calParams.pupil.X = meanP.X';
calParams.pupil.Y = meanP.Y';
calParams.glint.X = meanG.X';
calParams.glint.Y = meanG.Y';

%% load fixation points params
% load viewing distance
calParams.viewDist      = params.viewDist; % mm from screen

% load targets
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

%% Calculate the 'CalMat'
% Calculate the adjustment factor
calParams.rpc           = calcRpc(calParams);
% Calculate the transformation matrix
[calParams.calMat]      = calcCalMat(calParams);
calGaze                 = calcGaze(calParams);
theFig = figure;
hold on;

% plot each true and tracked target position. Red cross means target
% position and blue means tracked gaze position.
for i = 1:length(calGaze.X)
    plot(calParams.targets.X(i), calParams.targets.Y(i),'rx');
    plot(calGaze.X(i),calGaze.Y(i),'bx');
    plot([calParams.targets.X(i) calGaze.X(i)], [calParams.targets.Y(i) calGaze.Y(i)],'g');
end
errors = sqrt((calParams.targets.X-calGaze.X).^2+(calParams.targets.Y-calGaze.Y).^2);
accuracy = mean(errors(~isnan(errors)));
title(['Average error: ',num2str(accuracy),' mm'])
xlabel('Horizontal position (mm)');ylabel('Vertical position (mm)');
legend('Target position','Estimated gaze position')

% save the error figure
if ~exist(fullfile(dropboxDir, 'TOME_analysis', params.projectSubfolder, ...
        params.subjectName,params.sessionDate,'EyeTrackingQA'),'dir')
    mkdir(fullfile(dropboxDir, 'TOME_analysis', params.projectSubfolder, ...
        params.subjectName,params.sessionDate,'EyeTrackingQA'))
end
saveas(theFig,fullfile(dropboxDir, 'TOME_analysis', params.projectSubfolder, ...
    params.subjectName,params.sessionDate,'EyeTrackingQA', [params.calName '_MeanError']), 'pdf');
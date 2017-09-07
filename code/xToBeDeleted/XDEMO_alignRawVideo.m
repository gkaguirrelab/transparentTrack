%% Align raw tracking data to Report data

% This script aims to robustly align the RawTracking data to the TTL onset
% information provided in the LiveTrack Report.

% The LiveTrack system returns information about pupil size and position,
% gaze position (60 Hz) and TTL onset (30Hz) in a Report struct. The eye
% tracking data is not as accurate as we'd like to be, so we use a
% custom tracking routine to analyze the hi-res tracking video (raw video)
% acquired concurrently with the LiveTrack data. While the data
% obtained from the raw video (raw tracking) is much more accurate, it does
% not contain information about the TTL onset, and can't be directly synced
% to the TR timing.

% A property of the data acquisition system is that the raw video recording
% starts with a certain delay with respect to the Report recording. To
% align TTL information to the raw tracking data we need to calculate the
% data stream delay.

% We cannot rely on the timestamp provided by PsychToolbox for each field
% of the Report (i.e. frame) as an absolute time reference because it is
% assigned when the data is written in the variable, rather then when it is
% actually acquired by the device.

% We instead look at the signals resulting from the same measurement on the
% two data streams and make use of cross correlation (Signal Processing
% Toolbox) to align them and calculate the data stream delay.

%% Set defaults
% Get user name
[~, tmpName]            = system('whoami');
userName                = strtrim(tmpName);
% Set Dropbox directory
dbDir                   = ['/Users/' userName '/Dropbox-Aguirre-Brainard-Lab'];
% Set the subject / session / run
sessName                = 'session1_restAndStructure';
subjName                = 'TOME_3007';
sessDate                = '101116';
reportName              = 'rfMRI_REST_AP_run01_report.mat';
videoName               = 'rfMRI_REST_AP_run01_raw.mov';
outVideoFile            = fullfile('~','testVideo.avi');
outMatFile              = fullfile('~','testMat.mat');
numTRs                  = 420;
ltRes                   = [720 480] ./2; % rescaled resolution of the LiveTrack video (originally NTSC DV). Rescale is necessary because pupilTrack video is rescaled.
ptRes                   = [400 300]; % resolution of the pupilTrack video
ltThr                   = 0.1; % threshold for liveTrack glint position
ylims                   = [0.25 0.75];
%% Set the session and file names
sessDir                 = fullfile(dbDir,'TOME_data',sessName,subjName,sessDate,'EyeTracking');
reportFile              = fullfile(sessDir,reportName);
videoFile               = fullfile(sessDir,videoName);
%% Get the LiveTrack and raw video data
% LiveTrack
liveTrack               = load(reportFile);
% pupilTrack
params.inVideo          = videoFile;
params.outVideo         = outVideoFile;
params.outMat           = outMatFile;
[pupil,glint]           = trackPupil(params);
%% Perform a sanity checks on the LiveTrack report
% check if the frameCount is progressive
frameCountDiff          = unique(diff([liveTrack.Report.frameCount]));
assert(numel(frameCountDiff)==1,'LiveTrack frame Count is not progressive!');
% verify that the correct amount of TTLs has been recorded
[TTLPulses]             = CountTTLPulses (liveTrack.Report);
assert(TTLPulses==numTRs,'LiveTrack TTLs do not match TRs!');
% verify that the TTLs are correctly spaced, assuming that the acquisition
% rate is 30 Hz

%%% need to add this sanity check %%%

%% Use the X position of the glint to align data
% LiveTrack
%   average the two channels, output is at 30Hz
ltSignal                = mean([...
    [liveTrack.Report.Glint1CameraX_Ch01];...
    [liveTrack.Report.Glint1CameraX_Ch02]]);
ltNorm                  = ltSignal / ltRes(1);
% Remove poor tracking
ltDiff                  = [0 diff(ltNorm)];
ltNorm(abs(ltDiff) > ltThr)  = nan; % remove glint positions < ltThr
% pupilTrack
ptSignal                = glint.X;
ptNorm                  = (ptSignal / ptRes(1))';
%% Cross correlate the signals to compute the delay
% cross correlation doesn't work with NaNs, so we change them to zeros
ltCorr                  = ltNorm;
ptCorr                  = ptNorm;
ltCorr(isnan(ltNorm))   = 0 ;
ptCorr(isnan(ptNorm))   = 0 ;
% set vectors to be the same length
if length(ptCorr) > length(ltCorr)
    ltNorm              = [ltNorm,zeros(1,(length(ptCorr) - length(ltCorr)))];
    ltCorr              = [ltCorr,zeros(1,(length(ptCorr) - length(ltCorr)))];
else
    ptNorm              = [ptNorm,zeros(1,(length(ltCorr) - length(ptCorr)))];
    ptCorr              = [ptCorr,zeros(1,(length(ltCorr) - length(ptCorr)))];
end
% calculate cross correlation and lag array
[r,lag]                 = xcorr(ltCorr,ptCorr);

% when cross correlation of the signals is max the lag equals the delay
[~,I]                   = max(abs(r));
delay                   = lag(I); % unit = [number of samples]

% shift the signals by the 'delay'
ltAligned               = ltNorm;
ltAligned(ltAligned==0) = nan;
ptAligned               = [zeros(1,delay),ptNorm(1:end-delay)];
ptAligned(ptAligned==0) = nan;
%% Plot the results
fullFigure;
% before alignment
subplot(2,1,1)
plot(ltNorm, 'LineWidth',2);
hold on;
plot(ptNorm, 'LineWidth',2)
grid on
ylabel('glint X (normalized)')
xlabel('Frames')
legend ('liveTrack','pupilTrack')
title ('Before alignment')
ylim(ylims);
% after alignment
subplot(2,1,2);
plot(ltAligned, 'LineWidth',2);
hold on;
plot(ptAligned, 'LineWidth',2)
grid on
ylabel('glint X (normalized)')
xlabel('Frames')
legend ('liveTrack','pupilTrack')
title(['After alignment (shift = ' num2str(delay) ' frames)']);
ylim(ylims);
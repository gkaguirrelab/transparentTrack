function [delay] = calcRawVideoDelay(LiveTrackReport,RawTrackData, RawTrackRate)
% This function looks at the signals resulting from the same measurement
% made by the LiveTrack and by the custom tracking algorithm and makes use
% of cross correlation to align them and calculate the data stream delay.

% Usage
% 
% LiveTrackReport = load('/Users/giulia/Dropbox-Aguirre-Brainard-Lab/TOME_data/session1_restAndStructure/TOME_3004/091916/EyeTracking/rfMRI_REST_AP_run01_report.mat');
% RawTrackData = load('/Users/giulia/Dropbox-Aguirre-Brainard-Lab/eyeTrackingVideos/TOME_3004-rfMRI_REST_AP_run01.mat');
% RawTrackRate = 30; % in Hz
% [delay] = calcRawVideoDelay(LiveTrackReport,RawTrackData, RawTrackRate)


%% set default RawTrackRate
if ~exists ('RawTrackRate', 'var')
    RawTrackRate = 30;
end

%% Perform a sanity check on the LiveTrack Report
% check if the frameCount is progressive (i.e. no frame was skipped during
% data writing)
frameCountDiff = diff([LiveTrackReport.Report.frameCount]);
% note that there could be a frame count reset at the beginning, so it is
% ok if the frameCountDiff "jumps" once. If it jumps more, it is most
% likely an undesired skip.
frameCountJumps = diff(frameCountDiff);
skips = find(frameCountJumps);
if length(skips)>1
    warning('Frame Count is not progressive. Check current LiveTrack Report for integrity')
end

% further sanity checks for fMRI runs:
% verify that the correct amount of TTLs has been recorded
[TTLPulses] = CountTTLPulses (LiveTrackReport.Report);
fprintf(' >> LiveTrack recorded %d TTL pulses', TTLPulses);


%% Choose a reference signal to align data
% We need to choose the most precise LiveTrack measurement to help the
% aligning accuracy. For now, we use the X position of the glint by
% default.

% Raw Track signal
RTsignal = RawTrackData.glint.XY(1,:);

% Live Track signal
switch RawTrackRate
    case 30
        % if Raw tracking is done at 30 Hz, we average the data coming from
        % the two LiveTrack channels for each frame.
        LTsignal = ([LiveTrackReport.Report.Glint1CameraX_Ch01] + [LiveTrackReport.Report.Glint1CameraX_Ch02]) ./2;
    case 60
        % if Raw tracking is done at 60 Hz, we use all Report samples
        ct = 0;
        for i = 1:length(LiveTrackReport.Report)
            % First field
            ct                  = ct + 1;
            LTsignal(ct)         = LiveTrackReport.Report(i).Glint1CameraX_Ch01;
            % Second field
            ct                  = ct + 1;
            LTsignal(ct)         = LiveTrackReport.Report(i).Glint1CameraX_Ch02;
        end
end

%% Cross correlate the signals to compute the delay
% cross correlation doesn't work with NaNs, so we change them to zeros
RTsignal(isnan(RTsignal)) = 0 ;
LTsignal(isnan(LTsignal)) = 0 ;

% calculate cross correlation and also return the lag array
[r,lag] = xcorr(LTsignal,RTsignal);

% when cross correlation of the signals is max the lag equals the delay
[~,I] = max(abs(r));
delay = lag(I); % unit = [number of samples]


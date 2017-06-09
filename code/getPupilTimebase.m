function [timeBase] = getPupilTimebase(params,dropboxDir)

% This function assigns a timeBase to pupil data obtained using
% trackPupil.m and to liveTrack data. The timeBase is relative to the TTL
% pulsed recevied by the liveTrack device and is linked to the TR from the
% scanner and the 'T's on the stimulus computer. The liveTrack report and
% the pupil data are aligned cross-correlation the glint X position
% information on the 2 datasets (most reliable information from the
% LiveTrack report).
%
%   Usage:
%       [timebase] = getPupilTimebase(params, dropboxDir)
%
%   Required inputs:
%       params.outputDir
%       params.projectFolder
%       params.projectSubfolder
%       params.eyeTrackingDir
%
%       params.subjectName
%       params.sessionDate
%       params.runName
%       params.numTRs
%       params.acqRate
%       params.pupilTrackFile
%       params.ltRes
%       params.ptRes
%       params.maxLag
%       params.reportSanityCheck (default 1)
%       params.plotAlignment (default 0)
%
% Note that the params field are the same as the metaData fields for a
% standard pupilResponse struct, so this function can also be used like
% this:
%      [timeBase] = getPupilTimebase(metaData, dropboxDir)
%
%
%   Written by Andrew S Bock, Giulia Frazzetta - Nov.2016

%% DEMO - uncomment and run this session for an usage example
% %
% % Get user name
% [~, tmpName]            = system('whoami');
% userName                = strtrim(tmpName);
% % Set Dropbox directory
% dropboxDir                   = ['/Users/' userName '/Dropbox-Aguirre-Brainard-Lab'];
% % Set the subject / session / run
% params.outputDir        = 'TOME_processing';
% params.projectFolder    = 'TOME_processing';
% params.projectSubfolder = 'session1_restAndStructure';
% params.eyeTrackingDir   = 'EyeTracking';
% params.subjectName         = 'TOME_3011';
% params.sessionDate         = '111116';
% params.runName          = 'rfMRI_REST_AP_run01';
% params.numTRs           = 420;
% params.ltThr             = 0.1; % threshold for liveTrack glint position
% params.acqRate          = 60; % pupilTrack video frameRate in Hz
% params.pupilTrackFile   = ''; %%% <- NEED TO CHANGE THIS FOR DEMO MODE
% params.ltRes            = [360 240]; % resolution of the LiveTrack video (half original size)
% params.ptRes            = [400 300]; % resolution of the pupilTrack video


%% Set the session and file names
if isfield(params,'projectSubfolder')
    sessDir                 = fullfile(dropboxDir,params.projectFolder,...
        params.projectSubfolder,params.subjectName,params.sessionDate,params.eyeTrackingDir);
    reportFile              = fullfile(sessDir,[params.runName '_report.mat']);
    
    pupilTrackFile          = fullfile(dropboxDir,params.outputDir,...
        params.projectSubfolder,params.subjectName,params.sessionDate,params.eyeTrackingDir,[params.runName '_pupilTrack.mat']);
else
    sessDir                 = fullfile(dropboxDir,params.projectFolder,...
        params.subjectName,params.sessionDate,params.eyeTrackingDir);
    reportFile              = fullfile(sessDir,[params.runName '_report.mat']);
    
    pupilTrackFile          = fullfile(dropboxDir,params.outputDir,...
        params.subjectName,params.sessionDate,params.eyeTrackingDir,[params.runName '_pupilTrack.mat']);
end

%% set default max lag, sanity check, plot alignment
if ~isfield(params, 'maxLag')
    params.maxLag = 1000;
end
if ~isfield(params, 'reportSanityCheck')
    params.reportSanityCheck = 1;
end
if ~isfield(params, 'plotAlignment')
    params.plotAlignment = 0;
end
%% load the LiveTrack and raw track data
% check if data exists
% if ~isfield(params,'pupilTrackFile')
%     error('No raw tracking data for this run. Run trackPupil.m')
% end
% load matFiles
liveTrack               = load(reportFile);
pupilTrack              = load(pupilTrackFile);

%% Perform some sanity checks on the LiveTrack report
if params.reportSanityCheck
    % check if the frameCount is progressive
    frameCountDiff          = unique(diff([liveTrack.Report.frameCount]));
    assert(numel(frameCountDiff)==1,'LiveTrack frame Count is not progressive!');
    
    % verify that the correct amount of TTLs has been recorded (for fMRI runs)
    if isfield(params,'numTRs')
        [TTLPulses]             = CountTTLPulses (liveTrack.Report);
        assert(TTLPulses==params.numTRs,'LiveTrack TTLs do not match TRs!');
        % verify that the TTLs are correctly spaced, assuming that the acquisition
        % rate is 30 Hz
        
        %%% need to add this sanity check %%%
    end
else
    warning ('No sanity check is being performed on the LiveTrack report.')
end

%% Use the X position of the glint to align data
% LiveTrack
switch params.acqRate
    case 30
        % average the two channels, output is at 30Hz
        ltSignal                = mean([...
            [liveTrack.Report.Glint1CameraX_Ch01];...
            [liveTrack.Report.Glint1CameraX_Ch02]]);
    case 60
        % use all Report samples
        ct = 0;
        for i = 1:length(liveTrack.Report)
            % First field
            ct                  = ct + 1;
            ltSignal(ct)         = liveTrack.Report(i).Glint1CameraX_Ch01;
            % Second field
            ct                  = ct + 1;
            ltSignal(ct)         = liveTrack.Report(i).Glint1CameraX_Ch02;
        end
end
ltNorm                  = ltSignal / params.ltRes(1);
% Remove poor tracking
ltDiff                  = [0 diff(ltNorm)];
ltNorm(abs(ltDiff) > params.ltThr)  = nan; % remove glint positions < ltThr
% pupilTrack
ptSignal                = pupilTrack.glint.X;
ptNorm                  = (ptSignal / params.ptRes(1))';
%% Cross correlate the signals to compute the delay
% cross correlation doesn't work with NaNs, so we change them to zeros
ltCorr                  = ltNorm;
ptCorr                  = ptNorm;
ltCorr(isnan(ltNorm))   = 0 ;
ptCorr(isnan(ptNorm))   = 0 ;
% set vectors to be the same length (zero pad the END of the shorter one)
if length(ptCorr) > length(ltCorr)
    ltNorm              = [ltNorm,zeros(1,(length(ptCorr) - length(ltCorr)))];
    ltCorr              = [ltCorr,zeros(1,(length(ptCorr) - length(ltCorr)))];
else
    ptNorm              = [ptNorm,zeros(1,(length(ltCorr) - length(ptCorr)))];
    ptCorr              = [ptCorr,zeros(1,(length(ltCorr) - length(ptCorr)))];
end
% calculate cross correlation and lag array
[r,lag]                 = xcorr(ltCorr,ptCorr,params.maxLag);

% when cross correlation of the signals is max the lag equals the delay
[~,I]                   = max(abs(r));
delay                   = lag(I); % unit = [number of samples]

% shift the signals by the 'delay'
ltAligned               = ltNorm; % lt is not shifted
ltAligned(ltAligned==0) = nan;
ptAligned               = [zeros(1,delay),ptNorm(1:end-delay)];
ptAligned(ptAligned==0) = nan;

%% assign a common timeBase
% since ltSignal was not shifted and ptSignal is now aligned, we can assign
% a common timeBase
timeBaseTMP                = 1:length(ltAligned);
% get the first from liveTrack.Report
allTTLs                 = find([liveTrack.Report.Digital_IO1] == 1);
% if present, set the first TR to time zero
if ~isempty(allTTLs)
    firstTR           = allTTLs(1);
    timeBase.lt       = (timeBaseTMP - firstTR) * (1/params.acqRate); %liveTrack timeBase in [sec]
else
    timeBase.lt       = (timeBaseTMP - 1) * (1/params.acqRate);
end
timeBasePtTMP       = timeBase.lt + delay * (1/params.acqRate); %pupilTrack timeBase in [sec]

% make timebase.pt as long as ptSignal
if length(timeBasePtTMP)<length(ptSignal)
    %add missing values
    ptPadding = [timeBasePtTMP(end)+(1/params.acqRate): (1/params.acqRate) :timeBasePtTMP(end)+((length(ptSignal)-length(timeBasePtTMP)*(1/params.acqRate)))];
    timeBase.pt = [timeBasePtTMP ptPadding];
elseif length(timeBasePtTMP)>length(ptSignal)
    %trim timeBase.pt
    timeBase.pt = timeBasePtTMP(1:length(ptSignal));
else
    timeBase.pt = timeBasePtTMP;
end

%% Plot the cross correlation results
if params.plotAlignment
    fullFigure;
    % before alignment
    ylims                   = [0.25 0.75];
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
end
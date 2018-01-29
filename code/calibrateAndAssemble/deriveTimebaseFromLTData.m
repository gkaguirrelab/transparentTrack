function [timebase] = deriveTimebaseFromLTData(glintFileName,ltReportFileName,varargin)
% Derives a timebase for eyetracking data acquired with the LiveTrack+VTop setup
% 
% Description:
%   this function is to be used with the LiveTrack+VTop setup and assignes
%   a timebase to the raw data by aligning the tracked glint position
%   obtained with the livetrack (which includes the TR times information)
%   and the tracked glint position obtained with the custom tracking
%   algorithm. Note that this step is only necessary because raw data
%   acquired with the LiveTrack+VTop setup lacks syncing information with
%   the scanner TRs. 
%   As the X glint position is on average the best tracking results that
%   the livetrack algorithm can provide, we cross-correlate that timeseries
%   with our own tracked glint X position and determine the delay (in
%   frames) between the start of the datastream from the VTop device with
%   respect to the LiveTrack datastream.
% 
%   Note that if no TR information is found in the LiveTrack Report file
%   (i.e. for anatomical runs), the alignment with the scanner data is not
%   possible.
% 
% Input (required)
%   glintFileName       - full path to the matFile with the glint tracking
%                         results.
%   ltReportFileName    - full path to the livetrack generated "report file"
%
% Optional key/value pairs (analysis)
%   'timebaseFileName'  - full path to the file that will contain the
%                         timebase information.
%   'maxLag'            - max lag allowed between the two glint timeseries
%   'numTRs'            - number of expeted TR for this run. This is used
%                         only for sanity check purposes. In principle,
%                         just a single TTL syncing signal is enough to
%                         align a full run.
%   'rawVidFrameRate'   - framerate of the raw video in units of frames per
%                         second (fps).
%   'ltDataThreshold'   - threshold to clean up the livetrack signal before
%                         cross correlation.
%   'reportSanityCheck' - if true, a sanity check on the livetrack report
%                         will be performed before starting the alignment.
%                         The sanity check verifies that the livetrack
%                         framecount is progressive and that the registered
%                         TR correspond to the expected amount.
%   'plotAlignment'     - if set to true, a plot to verify the alignement
%                         will be generated.
%
% Optional key/value pairs (environment)
%  'tbSnapshot'         - This should contain the output of the
%                         tbDeploymentSnapshot performed upon the result
%                         of the tbUse command. This documents the state
%                         of the system at the time of analysis.
%  'timestamp'          - AUTOMATIC; The current time and date
%  'username'           - AUTOMATIC; The user
%  'hostname'           - AUTOMATIC; The host
% 
% Output
%	timebase            - structure with fields that contain the timebase
%                         information in milliseconds, and a meta field.
%
%% input parser

p = inputParser; p.KeepUnmatched = true;

% Required
p.addRequired('glintFileName',@ischar);
p.addRequired('ltReportFileName',@ischar);

% Optional analysis parameters
p.addParameter('timebaseFileName', '', @ischar);
p.addParameter('maxLag',500, @isnumeric);
p.addParameter('numTRs',420, @isnumerical);
p.addParameter('rawVidFrameRate',60, @isnumeric);
p.addParameter('ltDataThreshold',0.1, @isnumeric);
p.addParameter('reportSanityCheck',true, @islogical);
p.addParameter('plotAlignment',false, @islogical);

% Optional display and I/O parameters
p.addParameter('verbosity','none', @ischar);

% Environment parameters
p.addParameter('tbSnapshot',[],@(x)(isempty(x) | isstruct(x)));
p.addParameter('timestamp',char(datetime('now')),@ischar);
p.addParameter('username',char(java.lang.System.getProperty('user.name')),@ischar);
p.addParameter('hostname',char(java.net.InetAddress.getLocalHost.getHostName),@ischar);

% parse
p.parse(glintFileName,ltReportFileName, varargin{:})


%% load tracking data
glintData = load(glintFileName);
liveTrack = load(ltReportFileName);


%% Perform some sanity checks on the LiveTrack report
if p.Results.reportSanityCheck
    % check if the frameCount is progressive
    frameCountDiff = unique(diff([liveTrack.Report.frameCount]));
    if numel(frameCountDiff)==1
        warning('LiveTrack frame Count is not progressive!');
    end
    % verify that the correct amount of TTLs has been recorded (for fMRI runs)
    if ~isempty(p.Results.numTRs)
        allPulses = find ([liveTrack.Report.Digital_IO1]);
        if isempty(allPulses)
            TTLPulses = 0;
        else
            spacing = diff(allPulses);
            for ii = 1:length(spacing)
                if spacing (ii) == 1
                    adjacent (ii) = 1;
                else
                    adjacent (ii) = 0;
                end
            end
            TTLPulses = length (allPulses) - length (find(adjacent));
        end
        assert(TTLPulses==p.Results.numTRs,'LiveTrack TTLs do not match TRs!');
    end
else
    warning ('No sanity check is being performed on the LiveTrack report.')
end

%% Prepare data for alignment
% pick the right livetrack sampling data, based on the resolution of the
% raw video.
switch p.Results.rawVidFrameRate
    case 30
        % average the two channels, output is at 30Hz
        ltSignal = mean([[liveTrack.Report.Glint1CameraX_Ch01];...
            [liveTrack.Report.Glint1CameraX_Ch02]]);
    case 60
        % use all Report samples
        ct = 0;
        for ii = 1:length(liveTrack.Report)
            % First field
            ct = ct + 1;
            ltSignal(ct) = liveTrack.Report(ii).Glint1CameraX_Ch01;
            % Second field
            ct = ct + 1;
            ltSignal(ct) = liveTrack.Report(ii).Glint1CameraX_Ch02;
        end
end

% extract glint X position
glintSignal = glintData.glintData.X;


%% Cross correlate the signals to compute the delay
% cross correlation doesn't work with NaNs, so we change them to zeros
ltCorr = ltSignal;
glintCorr = glintSignal;
ltCorr(isnan(ltSignal)) = 0;
glintCorr(isnan(glintSignal)) = 0;

% set vectors to be the same length (zero pad the END of the shorter one)
if length(glintCorr) > length(ltCorr)
    ltSignal = [ltSignal,zeros(1,(length(glintCorr) - length(ltCorr)))];
    ltCorr = [ltCorr,zeros(1,(length(glintCorr) - length(ltCorr)))];
else
    glintSignal = [glintSignal; zeros((length(ltCorr) - length(glintCorr)),1)];
    glintCorr = [glintCorr; zeros((length(ltCorr) - length(glintCorr)),1)];
end

% calculate cross correlation and lag array
[r,lag]  = xcorr(ltCorr,glintCorr,p.Results.maxLag);

% when cross correlation of the signals is max the lag equals the delay
[~,I] = max(abs(r));
delay = lag(I); % unit = [number of samples]

% shift the signals by the 'delay' and replace the nans
ltAligned = ltSignal; % lt is not shifted
ltAligned(ltAligned==0) = nan;
glintAligned = [zeros(delay,1);glintSignal(1:end-delay)];
glintAligned(glintAligned==0) = nan;


%% assign a common timebase
% since ltSignal was not shifted and ptSignal is now aligned, we can assign
% a common timeBase
timebaseTMP = 1:length(ltAligned);

% get the TRs from liveTrack.Report
allTTLs = find([liveTrack.Report.Digital_IO1] == 1);

% if TR are present, set the first TR to time zero. Otherwise, set as zero
% the first liveTrack sample collected.
if ~isempty(allTTLs)
    firstTR = allTTLs(1);
    liveTrackTimebase = ((timebaseTMP - firstTR) * (1/p.Results.rawVidFrameRate))'; %liveTrack timeBase in [sec]
else
    liveTrackTimebase = ((timebaseTMP - 1) * (1/p.Results.rawVidFrameRate))';
end
timebaseGlintTMP = liveTrackTimebase + delay * (1/p.Results.rawVidFrameRate); %pupilTrack timeBase in [sec]

% make timebase.timebase as long as ptSignal
if length(timebaseGlintTMP)<length(glintData.glintData.X)
    %add missing values
    glintPadding = timebaseGlintTMP(end)+(1/p.Results.rawVidFrameRate): (1/p.Results.rawVidFrameRate) :timebaseGlintTMP(end)+((length(glintData.glintData.X)-length(timebaseGlintTMP)*(1/p.Results.rawVidFrameRate)));
    timebase.timebase = [timebaseGlintTMP glintPadding];
elseif length(timebaseGlintTMP)>length(glintData.glintData.X)
    %trim timeBase.rawVid
    timebase.timebase = timebaseGlintTMP(1:length(glintData.glintData.X));
else
    timebase.timebase = timebaseGlintTMP;
end

% convert in milliseconds
liveTrackTimebase = liveTrackTimebase * 1000;
timebase.timebase = timebase.timebase * 1000;
%% If so requested, plot the cross correlation results for quick review
if p.Results.plotAlignment
    figure;
    % before alignment
    subplot(2,1,1)
    plot(ltSignal, 'LineWidth',2);
    hold on;
    plot(glintSignal, 'LineWidth',2)
    grid on
    ylabel('glint X (different resolutions)')
    xlabel('Frames')
    legend ('liveTrack','glint tracked on raw video')
    title ('Before alignment')
    % after alignment
    subplot(2,1,2);
    plot(ltAligned, 'LineWidth',2);
    hold on;
    plot(glintAligned, 'LineWidth',2)
    grid on
    ylabel('glint X (different resolutions)')
    xlabel('Frames')
    legend ('liveTrack','glint tracked on raw video')
    title(['After alignment (shift = ' num2str(delay) ' frames)']);
end


%% save out the timebase data
% add some metafields first
timebase.meta = p.Results;
timebase.meta.delay = delay;
timebase.meta.units = 'milliseconds';
timebase.meta.liveTrackTimebase = liveTrackTimebase;

if ~isempty(p.Results.timebaseFileName)
    save(p.Results.timebaseFileName,'timebase');
end


end % main function

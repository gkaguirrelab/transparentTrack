function deriveTimebaseFromLTData(glintFileName,ltReportFileName,timebaseFileName,varargin)

% deriveTimebaseFromLTData(glintFileName,ltReportFilename,timebaseFilename

% header

%% input parser

p = inputParser; p.KeepUnmatched = true;

% Required
p.addRequired('glintFileName',@ischar);
p.addRequired('ltReportFileName',@ischar);
p.addRequired('timebaseFileName',@ischar);

% Optional analysis parameters
p.addParameter('maxLag',500, @isnumeric);
p.addParameter('reportSanityCheck',true, @islogical);
p.addParameter('plotAlignment',false, @islogical);
p.addParameter('numTRs',420, @isnumerical || @isempty);
p.addParameter('rawVidFrameRate',60, @isnumeric);
p.addParameter('ltDataThreshold',0.1, @isnumeric);

% Optional display and I/O parameters
p.addParameter('verbosity','none', @ischar);

% Environment parameters
p.addParameter('tbSnapshot',[],@(x)(isempty(x) | isstruct(x)));
p.addParameter('timestamp',char(datetime('now')),@ischar);
p.addParameter('username',char(java.lang.System.getProperty('user.name')),@ischar);
p.addParameter('hostname',char(java.net.InetAddress.getLocalHost.getHostName),@ischar);

% parse
p.parse(glintFileName,ltReportFileName,timebaseFileName, varargin{:})


%% load tracking data
glintData = load(glintFileName);
liveTrack = load(ltReportFileName);


%% Perform some sanity checks on the LiveTrack report
if p.Results.reportSanityCheck
    % check if the frameCount is progressive
    frameCountDiff = unique(diff([liveTrack.Report.frameCount]));
    assert(numel(frameCountDiff)==1,'LiveTrack frame Count is not progressive!');
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
        assert(TTLPulses==params.numTRs,'LiveTrack TTLs do not match TRs!');
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
glintSignal = glintData.glint.X;


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
    glintSignal = [glintSignal,zeros(1,(length(ltCorr) - length(glintCorr)))];
    glintCorr = [glintCorr,zeros(1,(length(ltCorr) - length(glintCorr)))];
end

% calculate cross correlation and lag array
[r,lag]  = xcorr(ltCorr,glintCorr,params.maxLag);

% when cross correlation of the signals is max the lag equals the delay
[~,I] = max(abs(r));
delay = lag(I); % unit = [number of samples]

% shift the signals by the 'delay' and replace the nans
ltAligned = ltSignal; % lt is not shifted
ltAligned(ltAligned==0) = nan;
glintAligned = [zeros(1,delay),glintSignal(1:end-delay)];
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
    timebase.lt = (timebaseTMP - firstTR) * (1/p.Results.rawVidFrameRate); %liveTrack timeBase in [sec]
else
    timebase.lt = (timebaseTMP - 1) * (1/p.Results.rawVidFrameRate);
end
timebaseGlintTMP = timebase.lt + delay * (1/p.Results.rawVidFrameRate); %pupilTrack timeBase in [sec]

% make timebase.rawVideo as long as ptSignal
if length(timebaseGlintTMP)<length(glintSignal)
    %add missing values
    glintPadding = timebaseGlintTMP(end)+(1/p.Results.rawVidFrameRate): (1/p.Results.rawVidFrameRate) :timebaseGlintTMP(end)+((length(glintSignal)-length(timebaseGlintTMP)*(1/p.Results.rawVidFrameRate)));
    timebase.rawVideo = [timebaseGlintTMP glintPadding];
elseif length(timebaseGlintTMP)>length(glintSignal)
    %trim timeBase.rawVid
    timebase.rawVideo = timebaseGlintTMP(1:length(glintSignal));
else
    timebase.rawVideo = timebaseGlintTMP;
end


%% If so requested, plot the cross correlation results for quick review
if p.Results.plotAlignment
    fullFigure;
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

save(timebaseFileName,'timebase');


end % main function

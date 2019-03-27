function makeTimebase(rawVideoInFileName, timebaseFileName, varargin)
% Derive timebase from raw video source, adjust for TTL audio signals
%
% Syntax:
%  makeTimebase(rawVideoInFileName, timebaseFileName, varargin)
% 
% Description:
%   This function finds the number of audio pulses converted from TTL 
%   signals and saves the timebase of those pulses.
%
% Inputs:
%   rawVideoInFileName    - Char vector. Full path to the raw pupil video
%   timebaseFileName      - Char vector. Full path to location where the 
%                           timebase is saved
%
% Optional key/value pairs (analysis):
%  'deinterlace'          - Logical. If set to true, the number of frames
%                           in the video is doubled to create the timebase.
%  'audioTrackSync'       - Logical. If set to true, the timing of audio
%                           spikes is used to synchronize the timebase.
%  'checkCountTRs'        - Scalar. Number of TRs 
%  'ttlAudioThreshScaler' - Scalar. Every value above this threshold is
%                           considered as TTL. If left to empty, the value
%                           will be set to the max signal in the audio
%                           channel divided by 2.5.
%  'proximityThresh'      - Scalar, in units of number of audio samples. If
%                           supra-threshold audio values are closer than
%                           this threshold, the multiple supra-threshold
%                           values are treated as a single TTL.
%  'makePlots'            - Logical, creates diagnostic plots. Default  
%                           false
%
% Optional key/value pairs (flow control):
%  'nFrames'              - Scalar. Analyze fewer than the total number of
%                           frames
%  'startFrame'           - Scalar. The frame on which to start.
%
% Examples:
%{
    videoIn = '/Users/aguirre/Dropbox (Aguirre-Brainard Lab)/MELA_data/Experiments/OLApproach_TrialSequenceMR/MRMMTMel/DataFiles/HERO_jxv1/2019-03-16/session_1/tfMRI_OLBLOCK_AP_Run1.mov'
    timbaseOut = '~/Desktop/tfMRI_OLBLOCK_AP_Run1_timebase.mat'
    makeTimebase( videoIn, timbaseOut)
%}

%% Parse vargin for options passed here
p = inputParser; p.KeepUnmatched = true;

% Required
p.addRequired('rawVideoInFileName',@isstr);
p.addRequired('timebaseFileName',@isstr);

% Optional parameters
p.addParameter('deinterlace',true,@islogical);
p.addParameter('audioTrackSync',false,@islogical);
p.addParameter('checkCountTRs',[], @(x)(isempty(x) | isscalar(x)));
p.addParameter('ttlAudioThreshScaler',0.4, @isscalar);
p.addParameter('proximityThresh',10,@isscalar);
p.addParameter('makePlots',false,@islogical);

% Optional flow control params
p.addParameter('nFrames',Inf,@isscalar);
p.addParameter('startFrame',1,@isscalar);

% parse
p.parse(rawVideoInFileName, timebaseFileName, varargin{:})


%% Get the frame timebase from raw video

% Prepare the video object
videoInObj = videoIOWrapper(rawVideoInFileName,'ioAction','read');

% Obtain the video parameters
if p.Results.nFrames == Inf
    nFrames = floor(videoInObj.Duration*videoInObj.FrameRate);
else
    nFrames = p.Results.nFrames;
end

% Clean up the horrible DropBox directory name
if ismac
    cleanFileName = rawVideoInFileName;
    cleanFileName = strrep(cleanFileName,' ','\ ');
    cleanFileName = strrep(cleanFileName,'(','\(');
    cleanFileName = strrep(cleanFileName,')','\)');
    sysCommand = ['GetFileInfo -d ' cleanFileName];
    [~,videoCreationDateTime] = system(sysCommand);
    % Remove the trailing carriage return
    videoCreationDateTime = videoCreationDateTime(1:end-1);
else
    videoCreationDateTime='';
end

% Obtain the frame rate and video frame duration
frameRate = videoInObj.FrameRate;
if p.Results.deinterlace
    frameRate = frameRate * 2;
end
videoFrameDur = (1/frameRate)*1000;

% Create the timebase values and meta data
timebase.values = ((p.Results.startFrame-1)*videoFrameDur:videoFrameDur:((p.Results.startFrame+nFrames-1)*2-1)*videoFrameDur)';
timebase.meta = p.Results;
timebase.meta.video.videonFrames = nFrames;
timebase.meta.video.videoframeDuration = videoFrameDur;
timebase.meta.video.videoFrameRate = frameRate;
timebase.meta.video.units = 'milliseconds';
timebase.meta.video.videoCreationDateTime = videoCreationDateTime;


%% Load audio channel
if p.Results.audioTrackSync
    % We store TTL timing information on the 1st channel of the audio file
    % associated with a video. A device is used to step down the TTL signal to
    % an audio input. This produces an audio "spike" in the channel data.
    [temp, audioFramesPerSec] = audioread(rawVideoInFileName);
    audioChannelData = temp(:,1);
    
    
    %% Calculate a default ttlAudioPulseThresh
    % The threshold is set by the maximum value found in the audio channel,
    % multipled by the ttlAudioThreshScaler
    ttlAudioPulseThresh = max(audioChannelData)*p.Results.ttlAudioThreshScaler;
    
    % Derive the frameDurationMsecs
    audioFrameDurationMsecs = 1/(audioFramesPerSec/1000);
    
    
    %% Find peaks
    [~, peakLocs] = findpeaks(audioChannelData,'MinPeakHeight',ttlAudioPulseThresh);
    
    % Remove any peak locations that are too close together
    peakLocs = peakLocs(~[diff(peakLocs)<p.Results.proximityThresh; false]);
    
    % If checkCountTRs is set, evaluate if the right number of TTL pulses was
    % received
    if ~isempty(p.Results.checkCountTRs)
        if length(peakLocs) ~= p.Results.checkCountTRs
            error('makeTimebase:failedCheckCountTRs','Observed number of TTL pulses is different than expected number of TRs. Make sure you entered the correct TR number or try changing the thresholds')
        end
    end
    
    % Obtain the time (in msecs) of the onset of the first TTL pulse
    firstTTLtimeMsecs = peakLocs(1) *audioFrameDurationMsecs;
    
    % Shift the timebase
    timebase.values = timebase.values - firstTTLtimeMsecs;
    
    % Update the timebase meta information
    timebase.meta.audio.audioFrameRate = audioFramesPerSec/1000;
    timebase.meta.audio.audioframeDuration = audioFrameDurationMsecs;
    timebase.meta.audio.units = 'milliseconds';

% Create diagnostic plot
if p.Results.makePlots
    figure
    x = 0:audioFrameDurationMsecs:(length(audioChannelData)-1)*audioFrameDurationMsecs;
    plot(x,audioChannelData,'-k')
    hold on
    plot(x(peakLocs),audioChannelData(peakLocs),'*r')
end

end

% Save the timebase file
save(p.Results.timebaseFileName,'timebase');

end


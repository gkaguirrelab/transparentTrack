function calcGazeData(pupilFileName,glintFileName,targetsFileName,gazeDataFileName,varargin)

% calcGazeData(pupilFileName,glintFileName,targetsFileName,gazeDataFileName)
%
% this function will calculate the data necessary to compute the
% gazeCalibration factors starting from pupil, glint and target files.

% TARGET ALIGNMENT
%
% If the presentation time of the target in each location was recorded
% along with the dot locations, the routine will align the raw data with
% the following strategy:
% 1.build a target timeseries in pseudo-pixel units: map the target
% positions as if they were position of an ideal pupil center on the frame.
% 2. cross correlate the pupil and target position timeseries.
% 3. cross correlate the pupil and target velocity timeseries to refine
% the results.

% If there is no information about the presentation time of the target, the
% routine will attempt to automatically align the pupil timeseries with the
% following strategy:
% 1. get the pseudo-velocity peaks from the target locations
% 2. extract velocity peaks from pupil timeseries
% 3. attempt cross correlation of the peaks
% 4. create a candidate target position timeseries
% 5. attempt cross correlation between pupil and target timeseries
% 6. if cross correlation is not good, try correcting the candidate target
% position timeseries, shifting to the next available pupil velocity peak.
% This strategy is not guaranteed to work, and it is recommended to review
% the alignment results before proceeding. In case the automatic alignment
% fails, it is usually possible to align the pupil timeseries by hand, as
% the fixations can usually be spotted by eye.
%
%
% OUTPUTS: (saved to file)
%   gazeCalData
%       targets = location of each target espressed as mm
%           on screen (origin center of the screen, X axis growing left to
%           right, Y axis growing top to bottom)
%       pupil = location of pupil center for each fixation in pixels (origin top left corner)
%       glint = location of glint center for each fixation in pixels
%       viewingDistance = viewing distance from the screen
%       meta
%           pupilTimeseries
%           targetPseudoTimeseries
%           pupilRaw
%           glintRaw
%           fixDurationSec
%           units
%
% INPUTS:
%   pupilFileName = name of the pupil file for the current calibration
%   glintFileName = name of the glint file for the current calibration
%   targetsFileName = name of the targets file for the current calibration.
%       must contain target location in screen coordinates, should contain
%       target onset times (the routine will attempt to align target and
%       pupil data with no onset times, but it is not guaranteed to work).
%   gazeDataFileName - name of the mat file to save the gaze calibration
%       data
%
%  Optional params:
%   viewingDistance - distance of subject's eye from stimulus screen
%   useLiveTrackGazeData - set to true to use LiveTrack data for calibration
%   rawDataPath: string to the raw data path, if different from the
%       location of the liveTrack data.
%   ltFileSuffixLength - number of character appended after the "gaze
%       calibration run" name (e.g. GazeCal4_LTdat.mat, suffix length = 10)
%   frameRate - frameRate of the raw video data.
%   medFilterOrder - order of the despiking filter applied to the pupil
%       position timeseries. Should be above 100 to get a clean signal.
%   saccadeDistance - estimate distance in frames between saccades. This is
%       used to refine velocity peak extraction
%   minDelay: minimum delay allowed for the alignement of pupil and target
%       timeseries
%  fixationWindowPct - percentage of the smallest fixation duration to
%  consider for the moving median window.
%
%
% Optional key/value pairs (display and I/O)
%  'verbosity' - level of verbosity. [none, full]
%
% Options (environment)
%   tbSnapshot - the passed tbSnapshot output that is to be saved along
%      with the data
%   timestamp / username / hostname - these are automatically derived and
%      saved within the p.Results structure.


%% input parser

p = inputParser; p.KeepUnmatched = true;

% Required
p.addRequired('pupilFileName',@ischar);
p.addRequired('glintFileName',@ischar);
p.addRequired('targetsFileName',@ischar);
p.addRequired('gazeDataFileName',@ischar);

% Optional analysis parameters
p.addParameter('units', 'mm', @ischar)
p.addParameter('dataIsAligned',false, @islogical);
p.addParameter('frameRate', 60, @isnumeric)
p.addParameter('medFilterOrder', 150, @isnumeric)
p.addParameter('saccadeDistance', 30, @isnumeric)
p.addParameter('minDelay', 450, @isnumeric)
p.addParameter('fixationWindowPct', 50, @isnumerc)

% Optional display and I/O parameters
p.addParameter('verbosity','none', @ischar);
p.addParameter('showFigures',false, @islogical);

% Environment parameters
p.addParameter('tbSnapshot',[],@(x)(isempty(x) | isstruct(x)));
p.addParameter('timestamp',char(datetime('now')),@ischar);
p.addParameter('username',char(java.lang.System.getProperty('user.name')),@ischar);
p.addParameter('hostname',char(java.net.InetAddress.getLocalHost.getHostName),@ischar);

% parse
p.parse(pupilFileName,glintFileName,targetsFileName,gazeDataFileName,varargin{:})


%% load data
% target
targetData = load(targetsFileName);
gazeCalData.targets.X = targetData.targets.X;
gazeCalData.targets.Y = targetData.targets.Y;
gazeCalData.targets.layout = targetData.targets.layout;
gazeCalData.viewingDistance = targetData.targets.viewingDistance;
if isfield(targetData.targets, 'onsetTimes')
    gazeCalData.targets.times = targetData.targets.times;
    gazeCalData.fixDurationSec = diff(targetData.targets.times);
    targetsTimesRecorded  = true;
else
    % we will need to estimate the duration of each fixation by the raw data
    targetsTimesRecorded = false;
end
clear targetData

% pupil
pupilData = load(pupilFileName);
pupil.X = pupilData.pupilData.pPosteriorMeanTransparent(:,1);
pupil.Y = pupilData.pupilData.pPosteriorMeanTransparent(:,2);
clear pupilData

% glint
glintData = load(glintFileName);
glint.X = glintData.glintData.X;
glint.Y = glintData.glintData.Y;
clear glintData


if ~p.Results.dataIsAligned
    %% convert target locations to pseudo-pixel coordinates
    % This conversion is necessary to cross correlate pupil and target signals
    % and align them. If the pupil and glint files are aligned in time with the
    % target presentation times, this might be skipped.
    switch gazeCalData.targets.layout
        case '3x3grid'
            % The "pseudo pixel" values are based on the assumption
            % that the central target would be "seen" on the frame in the central
            % position [180 120] (if res [320, 240]), and all other targets would be
            % displaced 20 pixels around it according to their location.
            
            % available targets locations
            highTRG = max(LTdata.targets(:,1));
            centerTRG = 0;
            lowTRG = min(LTdata.targets(:,1));
            
            % conversion table (note that the X must be flipped because the high res
            % video is acquired mirrored).
            centerPPX = [180 120];
            highPPX = [160 140];
            lowPPX = [200 100];
            
            targetPPX_X = gazeCalData.targets.X;
            targetPPX_X(targetPPX_X == centerTRG) = centerPPX(1);
            targetPPX_X(targetPPX_X == highTRG) = highPPX(1);
            targetPPX_X(targetPPX_X == lowTRG) = lowPPX(1);
            
            targetPPX_Y = gazeCalData.targets.Y;
            targetPPX_Y(targetPPX_Y == centerTRG) = centerPPX(2);
            targetPPX_Y(targetPPX_Y == highTRG) = highPPX(2);
            targetPPX_Y(targetPPX_Y == lowTRG) = lowPPX(2);
        otherwise
            error('This target layout does not currently support auto alignment')
    end
    
    % extract pupil position
    pupilPosition = sqrt((pupil.X).^2 + (pupil.Y).^2);
    
    % remove nans
    pupilPosition = fillmissing(pupilPosition,'linear');
    
    % filter with high-order to remove spikes from noises and saccades
    pupilPosition = medfilt1(pupilPosition,p.Results.medFilterOrder,'truncate');
    
    % derive pupil velocity
    pupilVelocity = diff(pupilPosition);
    
    
    % if the dot times were recorded, we build the target timeseries in
    % pseudo-pixels and do a 2 steps cross correlation using position and
    % velocity.
    if dotTimesRecorded
        % build targets timeseries
        target.X = [];
        target.Y = [];
        targetDurFrames = round(gazeCalData.fixDurationSec .* p.Results.frameRate) ;
        for cc = 1 : length(targetPPX_X)
            target.X = [target.X; targetPPX_X(cc) .* ones(targetDurFrames(cc),1)];
            target.Y = [target.Y; targetPPX_Y(cc) .* ones(targetDurFrames(cc),1)];
        end
        
        % extract target position and velocity
        targetPosition = sqrt((target.X).^2 + (target.Y).^2);
        targetVelocity = diff(sqrt((target.X).^2 + (target.Y).^2));
        
        % first pass: cross correlate positions with a lag window length of
        % half the target timeseries length.
        [r,lag]= xcorr(pupilPosition,targetPosition);%, round(length(targetPosition)/2));
        [~,I]= max(abs(r));
        delayP = lag(I);
        
        % trim the pupilPosition and the pupil velocity
        if delayP == 0
            pupilPositionTrim = pupilPosition;
        else
            pupilPositionTrim = pupilPosition(delayP:end);
        end
        pupilVelocityTrim = diff(pupilPositionTrim);
        
        % second pass: cross correlate the velocity profile to refine the
        % results
        [r,lag]= xcorr(pupilVelocityTrim,targetVelocity);
        [~,I]= max(abs(r));
        delayV = lag(I);
        
        % refine pupil and glint data
        delay = delayP + delayV;
        firstTargetOnsetIDX = delay;
        lastTargetOffsetIDX = delay+length(targetPosition);
        pupil.X = pupil.X(firstTargetOnsetIDX:lastTargetOffsetIDX);
        pupil.Y = pupil.Y(firstTargetOnsetIDX:lastTargetOffsetIDX);
        glint.X = glint.X(firstTargetOnsetIDX:lastTargetOffsetIDX);
        glint.Y = glint.Y(firstTargetOnsetIDX:lastTargetOffsetIDX);
        
    end
    
    
    % if the dot times were not recorded we try to match the velocity peaks of
    % the pupil timeseries with the pseudo velocity peaks of the moving
    % target. While we cannot build a proper timeseries for the moving target,
    % we know that the target moves only within a defined interval [minFixTime
    % maxFixTime]. This inteval is hardcoded in the gazeCalibration function.
    
    if ~dotTimesRecorded
        try
            % extract target velocity peaks and valleys
            targetVelocityPeaks = diff(sqrt((targetPPX_X).^2 + (targetPPX_Y).^2));
            
            % extract pupil velocity peaks and valleys indexes
            [~, vPeaksIDX] = findpeaks(pupilVelocity,'MinPeakHeight',std(pupilVelocity),'MinPeakDistance',p.Results.saccadeDistance);
            invertedPupilVelocity = - pupilVelocity;
            [~, vValleysIDX] = findpeaks(invertedPupilVelocity,'MinPeakHeight',std(pupilVelocity),'MinPeakDistance',p.Results.saccadeDistance);
            
            % build pupil velocity peak array
            pupilVelocityPeakIDX = sort([vPeaksIDX; vValleysIDX]);
            
            % find and remove big short saccades that survived filtering
            saccades = find (diff(pupilVelocityPeakIDX) < p.Results.saccadeDistance);
            
            % we define a saccade as 3 peaks too close to each other. In that
            % case, the "landing peak" (the last) is the one to be preserved
            if ~isempty(saccades)
                pupilVelocityPeakIDX(saccades) = [];
            end
            
            % get the true peaks
            pupilVelocityPeaks = pupilVelocity(pupilVelocityPeakIDX);
            
            % cross correlate the peak arrays: the delay index is onset of the
            % second target.
            [r,lag]= xcorr(pupilVelocityPeaks,targetVelocityPeaks);
            [~,I]= max(abs(r));
            delayV = lag(I);
            
            alignementCorrect = 0;
            
            while ~alignementCorrect
                
                firstTargetOnsetIDX = pupilVelocityPeakIDX(delayV-1);
                lastTargetOnsetIDX = pupilVelocityPeakIDX((delayV-1)+length(targetVelocityPeaks));
                
                % refine the tail of the pupilPosition, cropping it when the fixation
                % reasonably ends.
                pupilPositionTail = pupilPosition(lastTargetOnsetIDX:end);
                
                % get an estimate of the targets duration on screen based on the
                % peaks
                averageTargetDur = round(mean(diff(pupilVelocityPeakIDX((delayV-1):(delayV-1)+length(targetVelocityPeaks)))));
                
                % compute a moving std for the pupilPosition tail and select the most
                % stable window
                [~,minWindowIDX] = min(movstd(pupilPositionTail,averageTargetDur));
                
                tailIDX = round(minWindowIDX+(averageTargetDur/2));
                
                % get last target offset
                lastTargetOffsetIDX = min(lastTargetOnsetIDX+tailIDX,length(pupilPosition));
                
                % get target durations in frames
                if length(pupilVelocityPeakIDX) > (delayV)+length(targetVelocityPeaks)
                    targetDurFrames = diff(pupilVelocityPeakIDX((delayV-1):(delayV)+length(targetVelocityPeaks)));
                else
                    targetDurFrames = diff([pupilVelocityPeakIDX((delayV-1):(delayV-1)+length(targetVelocityPeaks)); length(pupilVelocity)]);
                end
                
                % build the candidate pseudo target timeseries
                target.X = [];
                target.Y = [];
                for cc = 1 : length(targetPPX_X)
                    target.X = [target.X; targetPPX_X(cc) .* ones(targetDurFrames(cc),1)];
                    target.Y = [target.Y; targetPPX_Y(cc) .* ones(targetDurFrames(cc),1)];
                end
                
                % extract target position
                targetPosition = sqrt((target.X).^2 + (target.Y).^2);
                
                % check if these correlate well. If not, the alignement failed.
                [r,lag]= xcorr(pupilPosition,targetPosition);
                [~,I]= max(abs(r));
                delayP = lag(I);
                
                if delayP <= 0 && length(pupilVelocityPeakIDX(delayV:end))  >= length(targetVelocityPeaks)
                    % remove first peak and shitft the delay to the next peak
                    delayV = delayV +1;
                elseif delayP <= p.Results.minDelay && length(pupilVelocityPeakIDX(delayV:end))  >= length(targetVelocityPeaks)
                    delayV = delayV +1;
                else
                    delay = delayP;
                    alignementCorrect = 1; % or stop trying
                end
            end
            
            % get the data
            firstTargetOnsetIDX = delay;
            lastTargetOffsetIDX = delay+length(targetPosition);
            pupil.X = pupil.X(firstTargetOnsetIDX:lastTargetOffsetIDX);
            pupil.Y = pupil.Y(firstTargetOnsetIDX:lastTargetOffsetIDX);
            glint.X = rawGlintData.glintData.X(firstTargetOnsetIDX:lastTargetOffsetIDX);
            glint.Y = rawGlintData.glintData.Y(firstTargetOnsetIDX:lastTargetOffsetIDX);
            
            % also, record the estimated target duration in seconds
            gazeCalData.fixDurationSec = targetDurFrames ./p.Results.frameRate;
            
        catch ME
            figure
            plot(pupilPosition)
            title('Raw Pupil position timeseries')
            xlabel('Frames')
            ylabel('Position')
            error (['Automatic alignement failed. Please check the calibration data manually.' ME])
        end
    end
    
    %% show the raw pupil position and the estimated fixation durations
    % to confirm they are correct
    if p.Results.showFigures
        figure
        plot(pupilPosition(firstTargetOnsetIDX:lastTargetOffsetIDX))
        hold on
        plot(targetPosition)
        legend('Pupil','Targets')
        title('Alignement of pupil center and targets timeseries')
        xlabel('Frames')
        ylabel('Position (pixels)')
    end
    
else
    % if all data is aligned, just pull the target duration in frames
     targetDurFrames = round(gazeCalData.fixDurationSec .* p.Results.frameRate) ;
end

%% Get mean pupil and glint position for each target fixation

% make index vector based on target duration
dataIDX = 1;
for ct = 1 :length(gazeCalData.targets.X)
    dataIDX(ct+1) = dataIDX(ct) + targetDurFrames(ct);
end

% make moving median on window with with a percentage of the frames for
% each fixation
windowDur = min(round(targetDurFrames ./100 .* p.Results.fixationWindowPct));

% compute moving medians and standard deviations
for ct = 1 :length(gazeCalData.targets.X)
    % medians
    movMedianP.X= movmedian(pupil.X(dataIDX(ct) :dataIDX(ct+1)), windowDur,'omitnan','Endpoints','discard');
    movMedianP.Y = movmedian(pupil.Y(dataIDX(ct) :dataIDX(ct+1)), windowDur,'omitnan','Endpoints','discard');
    movMedianG.X = movmedian(glint.X(dataIDX(ct) :dataIDX(ct+1)), windowDur,'omitnan','Endpoints','discard');
    movMedianG.Y = movmedian(glint.Y(dataIDX(ct) :dataIDX(ct+1)), windowDur,'omitnan','Endpoints','discard');
    % standard deviations
    movStdP.X = movstd(pupil.X(dataIDX(ct) :dataIDX(ct+1)), windowDur,'omitnan','Endpoints','discard');
    movStdP.Y = movstd(pupil.Y(dataIDX(ct) :dataIDX(ct+1)), windowDur,'omitnan','Endpoints','discard');
    movStdG.X = movstd(glint.X(dataIDX(ct) :dataIDX(ct+1)), windowDur,'omitnan','Endpoints','discard');
    movStdG.Y = movstd(glint.Y(dataIDX(ct) :dataIDX(ct+1)), windowDur,'omitnan','Endpoints','discard');
    % get minimum std for the pupil
    [minPstd,minPstdIDX] = min(abs(movStdP.X)+abs(movStdP.Y));
    % get the window that minimizes the std for each target
    gazeCalData.pupil.X(ct) = movMedianP.X(minPstdIDX);
    gazeCalData.pupil.Y(ct) = movMedianP.Y(minPstdIDX);
    gazeCalData.glint.X(ct) = movMedianG.X(minPstdIDX);
    gazeCalData.glint.Y(ct) = movMedianG.Y(minPstdIDX);
end


%%  add meta fields and save results
gazeCalData.meta = p.Results;
gazeCalData.meta.pupilTimeseries = pupil;
gazeCalData.meta.glintTimeseries = glint;
gazeCalData.meta.targetPseudoTimeseries = target;
gazeCalData.meta.pupilRaw.X = rawPupilData.pupilData.pPosteriorMeanTransparent(:,1);
gazeCalData.meta.pupilRaw.Y = rawPupilData.pupilData.pPosteriorMeanTransparent(:,2);
gazeCalData.meta.glintRaw = rawGlintData.glintData;

% save results
save (gazeDataFileName , 'gazeCalData')

end % main function

    

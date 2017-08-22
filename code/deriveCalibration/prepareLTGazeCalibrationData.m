function prepareLTGazeCalibrationData (LTdatFileName,gazeDataFileName,varargin)

% prepareGazeCalibrationData (LTdatFileName,gazeDataFileName)

% this function is to use for gaze calibration data collected with the
% LiveTrack device alone or with the Livetrack+VTop eye tracking system. It
% will export a gazeData file containing all information needed to
% calculate the gaze calibration conversion factors.
% 
% gazeCalData 
%       targets.X = location on X axes with origin in the center of the 
%           screen and growing left to right espressed as mm on screen
%       targets.Y = location on Y axes with origin in the center of the 
%           screen and growing top to bottom espressed as mm on screen
%       pupil.X = raw location of pupil center in pixels on the frame (origin top left corner)  
%       pupil.Y = raw location of pupil center in pixels on the frame (origin top left corner)
%       glint.X = raw location of glint center in pixels on the frame (origin top left corner)
%       glint.Y = raw location of glint center in pixels on the frame (origin top left corner)
%       viewingDistance = viewing distance from the screen in mm
%       meta
%           pupilRaw
%           glintRaw
%           fixDurationSec
%           


%% input parser

p = inputParser; p.KeepUnmatched = true;

% Required
p.addRequired('LTdatFileName',@ischar);
p.addRequired('gazeDataFileName',@ischar);

% Optional analysis parameters
p.addParameter('useLiveTrackGazeData',false, @islogical)
p.addParameter('ltFileSuffixLength', 10, @isnumeric)
p.addParameter('frameRate', 60, @isnumeric)
p.addParameter('movingMeanWindow', 1/3, @isnumerc)
p.addParameter('removeBoundaryFrames', 1/10, @isnumeric)
p.addParameter('viewingDistance', 1065, @isnumeric)

% Optional display and I/O parameters
p.addParameter('verbosity','none', @ischar);

% Environment parameters
p.addParameter('tbSnapshot',[],@(x)(isempty(x) | isstruct(x)));
p.addParameter('timestamp',char(datetime('now')),@ischar);
p.addParameter('username',char(java.lang.System.getProperty('user.name')),@ischar);
p.addParameter('hostname',char(java.net.InetAddress.getLocalHost.getHostName),@ischar);

% parse
p.parse(sizeDataFilesNames, sizeFactorsFileName, varargin{:})


%% load LTdat file
LTdata = load(LTdatFileName);


%% read data from ltDatFile and store them in gazeCalData

% 1. get target location
gazeCalData.targets.X     = LTdata.targets(:,1); % mm on screen, screen center = 0
gazeCalData.targets.Y     = LTdata.targets(:,2); % mm on screen, screen center = 0

% find and replace nan target, if any

% if the livetrack failed to track one target, it will be a NaN. Replace
% the NaN with the missing target location. NOTE THAT THIS ASSUMES THAT
% ONLY ONE TARGET IS MISSING

nanIDX = find(isnan(gazeCalData.targets.X));

if ~isempty(nanIDX)
   if length(nanIDX) > 1
       error ('There is more than one NaN target! Calibration might fail.')
   end
    % available targets locations
    highTRG = max(LTdata.targets(:,1));
    centerTRG = 0;
    lowTRG = min(LTdata.targets(:,1));
    
    allLocations = [...
        highTRG highTRG; ...
        highTRG centerTRG; ...
        highTRG lowTRG; ...
        centerTRG highTRG; ...
        centerTRG centerTRG; ...
        centerTRG lowTRG; ...
        lowTRG highTRG; ...
        lowTRG centerTRG; ...
        lowTRG lowTRG; ...
        ];
    
    % find value of the nan target (is the one missing when comparing the
    % targets to allLocations)
    missingTarget = find(~ismember(allLocations,[gazeCalData.targets.X gazeCalData.targets.Y], 'rows'));
    
    % replace NaN target value with missing value
    gazeCalData.targets.X(nanIDX) = allLocations(missingTarget, 1);
    gazeCalData.targets.Y(nanIDX) = allLocations(missingTarget, 2);
end

% get viewing distance
gazeCalData.viewingDistance = p.Results.viewingDistance;


%%  If so requested, get pupil and glint coordinates as tracked by the LiveTrack and exit.
if p.Results.useLiveTrackGazeData
    gazeCalData.pupil.X = LTdata.pupil(:,1);
    gazeCalData.pupil.Y = LTdata.pupil(:,2);
    gazeCalData.glint.X = LTdata.glint(:,1);
    gazeCalData.glint.Y = LTdata.glint(:,2);
    
    % add metadata
    gazeCalData.meta = p.Results;
    
    % save results 
    save (gazeDataFileName , 'gazeCalData')
    
    % no other further step is necessary
    return
end


%% Load raw data and dot times

% get dot times (if available)
if isfield(LTdata,'dotTimes')
    gazeCalData.fixDurationSec = diff(LTdata.dotTimes);
    dotTimesRecorded  = true;
else
    % we will need to estimate the duration of each fixation by the raw data
    dotTimesRecorded = false;
end

% get calibration file root name
calFileRoot = ltDatFileName(1:end-p.Results.ltFileSuffixLength);

% get raw pupil and glint data
rawPupilData = load ([calFileRoot '_pupil.mat']);
rawGlintData = load ([calFileRoot '_glint.mat']);

% extract X and Y pupil timeseries from the posterior
pupil.X = rawPupilData.pupilData.pPosteriorMeanTransparent(:,1);
pupil.Y = rawPupilData.pupilData.pPosteriorMeanTransparent(:,2);


%% convert target locations to pseudo-pixel coordinates
% This conversion is necessary to cross correlate pupil and target signals
% and align them. The "pseudo pixel" values are based on the assumption
% that the central target would be "seen" on the frame in the central
% position [180 120] (if res [320, 240]), and all other targets would be
% displaced 20 pixels around it according to their location.

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


%% align target and pupil timeseries

% extract pupil position and velocity
pupilPosition = sqrt((pupil.X).^2 + (pupil.Y).^2);
pupilVelocity = diff(sqrt((pupil.X).^2 + (pupil.Y).^2));

% remove nan from pupil position and velocity
pupilPosition = fillmissing(pupilPosition,'linear');
pupilVelocity = fillmissing(pupilVelocity,'linear');

% if the dot times were recorded, we build the target timeseries in
% pseudo-pixels and do a 2 steps cross correlation using position and
% velocity.
if dotTimesRecorded
    % build targets timeseries
    target.X = [];
    target.Y = [];
    for cc = 1 : length(targetPPX_X)
        target.X = [target.X; targetPPX_X(cc) .* ones(round(gazeCalData.fixDurationSec(cc) * p.Results.frameRate),1)];
        target.Y = [target.Y; targetPPX_Y(cc) .* ones(round(gazeCalData.fixDurationSec(cc) * p.Results.frameRate),1)];
    end
    
    % extract target position and velocity
    targetPosition = sqrt((target.X).^2 + (target.Y).^2);
    targetVelocity = diff(sqrt((target.X).^2 + (target.Y).^2));
    
    % first pass: cross correlate positions with a lag window length of
    % half the target timeseries length. 
    [r,lag]= xcorr(pupilPosition,targetPosition, round(length(targetPosition)/2));
    [~,I]= max(abs(r));
    delayP = lag(I);
    
    % trim the pupilPosition and the pupil velocity
    pupilPositionTrim = pupilPosition(delayP:delayP+length(targetPosition));
    pupilVelocityTrim = diff(pupilPositionTrim);
    
    % second pass: cross correlate the velocity profile to refine the
    % results
    [r,lag]= xcorr(pupilVelocityTrim,targetVelocity);
    [~,I]= max(abs(r));
    delayV = lag(I);
    
    % refine pupil data
    delay = delayP + delayV;
    pupil.X = pupil.X(delay:delay+length(target.X));
    pupil.Y = pupil.Y(delay:delay+length(target.Y));
end


% if the dot times were not recorded...























%%  add meta fields and save results
gazeCalData.meta = p.Results;
gazeCal.meta.rawPupil = rawPupilData;
gazeCal.meta.rawGlint = rawGlintData;









%% FROM OLD FUNCTION


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


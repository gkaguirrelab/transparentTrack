function prepareGazeCalibrationData (LTdatFileName,gazeDataFileName,varargin)

% prepareGazeCalibrationData

% gazeCalData 
%       targets.X = location on X axes growing left to right espressed as mm on screen
%       targets.Y = location on Y axes growing top to bottom espressed as mm on screen
%       pupil.X = raw location of pupil center in pixels on the frame (origin top left corner)  
%       pupil.Y = raw location of pupil center in pixels on the frame (origin top left corner)
%       glint.X = raw location of glint center in pixels on the frame (origin top left corner)
%       glint.Y = raw location of glint center in pixels on the frame (origin top left corner)
%       viewingDistance = viewing distance from the screen in mm
%       meta
%           pupilRaw
%           glintRaw
%           targetOnsetTimes
%           


%% input parser

p = inputParser; p.KeepUnmatched = true;

% Required
p.addRequired('LTdatFileName',@ischar);
p.addRequired('gazeDataFileName',@ischar);

% Optional analysis parameters
p.addParameter('useLiveTrackGazeData',false, @islogical)
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
    
    missingTarget = find(~ismember(allLocations,[gazeCalData.targets.X gazeCalData.targets.Y], 'rows'));
    
    % replace NaN target value with missing value
    
    gazeCalData.targets.X(nanIDX) = allLocations(missingTarget, 1);
    gazeCalData.targets.Y(nanIDX) = allLocations(missingTarget, 2);
end

% get viewing distance
gazeCalData.viewingDistance = p.Results.viewingDistace;

% get pupil and glint coordinates as tracked by the LiveTrack
if p.Results.useLiveTrackGazeData
    gazeCalData.pupil.X = LTdata.pupil(:,1);
    gazeCalData.pupil.Y = LTdata.pupil(:,2);
    gazeCalData.glint.X = LTdata.glint(:,1);
    gazeCalData.glint.Y = LTdata.glint(:,2);
    
else
    % get dot times (if it was saved)
    
    % get raw pupil and glint data
    
    % align targets and raw data
    
    % get mean pupil and glint location
    
end

% add meta fields

gazeCalData.meta = p.Results;
if ~p.Results.useLiveTrackGazeData
%     gazeCal.meta.rawPupilX = 
%     gazeCal.meta.rawPupilY =
%     gazeCal.meta.rawGlintX = 
%     gazeCal.meta.rawGlintY =
end







%% FROM OLD FUNCTION
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


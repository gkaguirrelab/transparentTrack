function copyLTdatToGazeData(LTdatFileName,gazeDataFileName,varargin)
% copyLTdatToGazeData(LTdatFileName,gazeDataFileName)

% this function is to use for gaze calibration data collected with the
% LiveTrack device. If the user wishes to use the livetrack data for
% calibration, this routine export the data in a standard format for later
% calibration steps.
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
%           fixDurationSec
%           units
% 
% INPUTS:
%   LTdatFileName = name of the LiveTrack data file for the current
%       calibration
%   gazeDataFileName = name of the mat file to save the gaze calibration
%       data
% 
%  Optional params:
%   viewingDistance - distance of subject's eye from stimulus screen
%   units
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
p.addRequired('LTdatFileName',@ischar);
p.addRequired('gazeDataFileName',@ischar);

% Optional analysis parameters
p.addParameter('viewingDistanceMm', 1065, @isnumeric)
p.addParameter('units', 'mm', @ischar)

% Optional display and I/O parameters
p.addParameter('verbosity','none', @ischar);
p.addParameter('showFigures',false, @islogical);

% Environment parameters
p.addParameter('tbSnapshot',[],@(x)(isempty(x) | isstruct(x)));
p.addParameter('timestamp',char(datetime('now')),@ischar);
p.addParameter('username',char(java.lang.System.getProperty('user.name')),@ischar);
p.addParameter('hostname',char(java.net.InetAddress.getLocalHost.getHostName),@ischar);

% parse
p.parse(LTdatFileName, gazeDataFileName, varargin{:})


%% load LTdat file
LTdata = load(LTdatFileName);


%% read data from ltDatFile and store them in gazeCalData

% get target location
gazeCalData.targets.X     = LTdata.targets(:,1); % mm on screen, screen center = 0
gazeCalData.targets.Y     = LTdata.targets(:,2); % mm on screen, screen center = 0


% get pupil and glint location for each fixation
gazeCalData.pupil.X = LTdata.pupil(:,1);
gazeCalData.pupil.Y = LTdata.pupil(:,2);
gazeCalData.glint.X = LTdata.glint(:,1);
gazeCalData.glint.Y = LTdata.glint(:,2);

% add viewing distance
gazeCalData.viewingDistanceMm = p.Results.viewingDistanceMm;

% add metadata
gazeCalData.meta = p.Results;

% add an additional field to indicate that we are using livetrack data
gazeCalData.meta.usingLiveTrackData = true;

% save results
save (gazeDataFileName , 'gazeCalData')

end % main function

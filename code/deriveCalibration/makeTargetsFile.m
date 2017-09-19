function makeTargetsFile(targetsInfoFile,targetsFileName,varargin)

% makeTargetsFile(targetInfoFile,targetsFileName)
% 
%  This routine will extract the information about the location and timing
%  of the target during the fixation task necessary for gaze calibration.
% 
% OUTPUTS: (saved to file)
%   targets - struct containing the following target info
%         X
%         Y
%         sysClockSecsOnsets
%         sysClockSecsOffsets
%         viewingDistanceMm
%         meta
% 
% INPUTS:
%   targetsInfoFile - name of the file containing targets information
%   targetsFileName - name of the targets file for saving the output.
% 
%  Optional params:
%   targetsInfoFileType - type of file from which to pull the target info
%   targetsLayout: description of the target layout
%   viewingDistance - frameRate of the raw video data.
%   targetsUnits - units in which the target position is expressed.
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
p.addParameter('targetsInfoFileType','LiveTrack', @ischar) % alternative '3secTarget'
p.addParameter('viewingDistanceMm', 1065, @isnumeric)
p.addParameter('targetsPositionUnits','mmOnScreen',@ischar);

% Optional display and I/O parameters
p.addParameter('verbosity','none', @ischar);

% Environment parameters
p.addParameter('tbSnapshot',[],@(x)(isempty(x) | isstruct(x)));
p.addParameter('timestamp',char(datetime('now')),@ischar);
p.addParameter('username',char(java.lang.System.getProperty('user.name')),@ischar);
p.addParameter('hostname',char(java.net.InetAddress.getLocalHost.getHostName),@ischar);

% parse
p.parse(targetsInfoFile,targetsFileName,varargin{:})


%% load the target info file
targetsInfo = load(targetsInfoFile);


%% pull the target data according to the type of target info file
switch p.Results.targetsInfoFileType
    case 'LiveTrack'
        targets.X = targetsInfo.targets(:,1); % mm on screen, screen center = 0
        targets.Y = targetsInfo.targets(:,2); % mm on screen, screen center = 0
        
        % if the livetrack failed to track one target, it will be a NaN. Replace
        % the NaN with the missing target location. NOTE THAT THIS ASSUMES THAT
        % ONLY ONE TARGET IS MISSING
        
        nanIDX = find(isnan(targets.X));
        if ~isempty(nanIDX)
            if length(nanIDX) > 1
                error ('There is more than one NaN target! Calibration might fail.')
            end
            % available targets locations
            highTRG = max(targetsInfo.targets(:,1));
            centerTRG = 0;
            lowTRG = min(targetsInfo.targets(:,1));
            
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
            missingTarget = find(~ismember(allLocations,[targets.X targets.Y], 'rows'));
            
            % replace NaN target value with missing value
            targets.X(nanIDX) = allLocations(missingTarget, 1);
            targets.Y(nanIDX) = allLocations(missingTarget, 2);
        end
        
        % get dot times (if available)
        if isfield(targetsInfo,'dotTimes')
            targets.sysClockSecsOnsets = targetsInfo.dotTimes(1:end-1)';
            targets.sysClockSecsOffsets = targetsInfo.dotTimes(2:end)';
        end
        
        
    case '3secTarget'
        % get targets location
        targets.X     = targetsInfo.targets(:,1); % mm on screen, screen center = 0
        targets.Y     = targetsInfo .targets(:,2); % mm on screen, screen center = 0
           
        % get targets times
        targets.sysClockSecsOnsets = targetsInfo.dotTimes(1:end-1)';
        targets.sysClockSecsOffsets = targetsInfo.dotTimes(2:end)';
        
    otherwise
        error('Unknown targetsInfoFileType')
end


% get viewing distance and layout
targets.viewingDistanceMm = p.Results.viewingDistance;

%% add a meta field and save out the results

targets.meta = p.Results;

save(targetsFileName,'targets');


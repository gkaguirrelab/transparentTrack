function runCalibrationAndTimebasePipeline(pathParams, varargin )

% header

% steps:
% calcSizeCalFactors
% applySizeCalibration
% prepareLTGazeCalibrationData
% calcGazeCalFactors
% applyGazeCalibration
% deriveTimebaseFromLTData
%% Parse input and define variables
p = inputParser; p.KeepUnmatched = true;

% required input
p.addRequired('pathParams',@isstruct);

% optional input
p.addParameter('videoTypeChoice', 'LiveTrackWithVTOP', @ischar);
p.addParameter('variableNamingConvention', 'LiveTrackWithVTOP', @ischar);
p.addParameter('customFunCalls', {}, @iscell);
p.addParameter('skipStage', {}, @iscell);
p.addParameter('sizeCalIdentifier', '*Scale*_pupil.mat', @ischar);
p.addParameter('sizeCalSuffixLength', '10', @isnumeric);
p.addParameter('mostRecentGazeCal', 'before', @ischar); % alternative 'after'

% parse
p.parse(pathParams, varargin{:})
pathParams=p.Results.pathParams;

% Sanity check the input
if ~isempty(p.Results.customFunCalls) && ~strcmp(p.Results.videoTypeChoice,'custom')
    error('Set videoTypeChoice to ''custom'' when passing customFunCalls');
end

if strcmp(p.Results.videoTypeChoice,'custom') && isempty(p.Results.customFunCalls)
    error('customFunCalls is set to empty, despite request for a custom videoTypeChoice');
end

%% Create output directories if needed
if ~exist(pathParams.dataOutputDirFull,'dir')
    mkdir(pathParams.dataOutputDirFull)
end
if ~exist(pathParams.controlFileDirFull,'dir')
    mkdir(pathParams.controlFileDirFull)
end


%% Define input and output filenames
glintFileName = fullfile(pathParams.dataOutputDirFull, [pathParams.runName '_glint.mat']);
pupilFileName = fullfile(pathParams.dataOutputDirFull, [pathParams.runName '_pupil.mat']);

% the calibration routine might greatily vary with the acquisition method.
% We then use the variableNamingConveition as a switch to define the variable names.
switch p.Results.variableNamingConvention
    case 'LiveTrackWithVTOP'
        % for size factors
        allSizeRawFiles = dir(fullfile(pathParams.dataOutputDirFull, p.Results.sizeCalIdentifier));
        for ii = 1:length(allSizeRawFiles)
            sizeDataFilesNames{ii} = fullfile(pathParams.dataOutputDirFull,allSizeRawFiles(ii).name(1:end-p.Results.sizeCalSuffixLength));
        end
        
        % for size calibration
        sizeCalFactorsFileName = fullfile(pathParams.dataOutputDirFull, [pathParams.runName '_sizeCalFactors.mat']);
        calibratedPupilFileName = fullfile(pathParams.dataOutputDirFull, [pathParams.runName '_calibratedPupil.mat']);
        
        % for gaze factors
        LTdatFileName = pickLTGazeData(pathParams,p.Results.mostRecentGazeCal);
        gazeDataFileName = fullfile(pathParams.dataOutputDirFull, [pathParams.runName '_gazeCalData.mat']);
        gazeCalFactorsFileName = fullfile(pathParams.dataOutputDirFull, [pathParams.runName '_gazeCalFactors.mat']);
        
        % for gaze calibration
        gazeCalFactorsFileName = fullfile(pathParams.dataOutputDirFull, [pathParams.runName '_calibratedGaze.mat']);
        
        % for timebase
        ltReportFileName = fullfile(pathParams.dataSourceDirFull, [pathParams.runName '_report.mat']);
        timebaseFileName = fullfile(pathParams.dataOutputDirFull, [pathParams.runName '_timebase.mat']);
        
    otherwise
        error('Variable names not defined for this variableNamingConvention')
end

%% define analysis pipelines.
switch p.Results.videoTypeChoice
    case 'LiveTrackWithVTOP'
        funCalls = {...
            'deriveTimebaseFromLTData(glintFileName,ltReportFileName,timebaseFileName, varargin{:});'...
            'calcSizeCalFactors(sizeDataFilesNames, sizeCalFactorsFileName, varargin{:});'...
            'applySizeCalibration(pupilFileName,sizeFactorsFileName,calibratedPupilFileName, varargin{:});'...
            'prepareLTGazeCalibrationData (LTdatFileName,gazeDataFileName, varargin{:});'...
            'calcGazeCalFactors (gazeDataFileName,gazeCalFactorsFileName, varargin{:});' ...
            'applyGazeCalibration(pupilFileName,glintFileName,gazeCalFactorsFileName,calibratedGazeFileName, varargin{:});' ...
            };
    case 'custom'
        funCalls = p.Results.customFunCalls;
    otherwise
        error('VideoTypeChoice is not defined')
end

% Grab the function names as the portion of the funCalls that preceed the
% open parenthesis
funNames = cellfun(@(x) strtok(x,'('),funCalls,'UniformOutput',false);

% Loop through the function calls
for ff = 1:length(funCalls)
    if ~any(strcmp(p.Results.skipStage,funNames{ff}))
        eval(funCalls{ff});
        if strcmp(p.Results.lastStage,funNames{ff})
            return
        end
    end
end

end % main function


function [LTdatFileName] = pickLTGazeData(pathParams,mostRecentGazeCal)
% this function will automatically select the calibration that is most
% appropriate for the run

GazeData = dir(pathParams.dataSourceDirFull,'*LTdat*.mat');
if length(GazeCals) == 1
    LTdatFileName = fullfile(pathParams.dataSourceDirFull,GazeData.name);
elseif length(GazeCals) > 1
    % sort Gaze Data files by timestamp
    [~,idx] = sort([GazeData.datenum]);
    % get the timestamp for the current run from the corresponding
    % report
    reportFile = dir(fullfile(pathParams.dataSourceDirFull,[params.runName '_report.mat']));
    reportTime = reportFile.datenum;
    switch mostRecentGazeCal
        case 'before'
            % take the most recent calibration file acquired before the current run
            for ii = 1: length(idx)
                if GazeCals(idx(ii)).datenum < reportTime && isGood(idx(ii))
                    LTdatFileName = (fullfile(pathParams.dataSourceDirFull,GazeCals(idx(ii)).name));
                else
                    continue
                end
            end
    end
else
    warning('No gaze calibration found for this run. Skipping gaze calibration step')
    LTdatFileName = [];
end
end

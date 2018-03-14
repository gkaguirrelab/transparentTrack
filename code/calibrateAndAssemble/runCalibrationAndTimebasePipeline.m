function runCalibrationAndTimebasePipeline(pathParams, varargin )
% runs the full calibration pipeline
% 
% Description
%   This is a wrapper to run the full standard calibration pipeline for
%   videos acquired with the LiveTrack+VTop setup. It includes the
%   following steps:
%       deriveTimebaseFromLTData
%       calcSizeCalFactors
%       applySizeCalibration
%       makeTargetsFile
%       calcGazeCalFactors
%       applyGazeCalibration
%   Each step is fully documented in the header of the corresponding
%   function.
% 
% Input (required)
%   pathParams 			       - Structure with the path information of the 
%                                files to be analyzed.
% 
% Optional key/value pairs (analysis)
%   'videoTypeChoice'          - Type of video acquisition setup (only 
%                                option currently available is
%                                'LiveTrackWithVTOP', or the flag 'custom'
%                                can be used followed in combination with
%                                the 'customFunCall' option).
%   'variableNamingConvention' - Type of eytetracking file naming convention
%                                for automatic retreival (only
%                                optioncurrently available is
%                               'LiveTrackWithVTOP').
%   'customFunCalls'           - Cell array to customize the calls 
%                                to the calibration functions.
%   'skipStageByName'          - A cell array of function calls to be
%                                skipped during execution of the pipeline.
%   'lastStage'                - The last stage to be executed. By deafult
%                                ends with the application of the gaze
%                                calibration.
%   'sizeCalIdentifier'        - Combination of wildcart and characters
%                                used to automatically identify size
%                                calibration files.
%   'sizeCalSuffixLength'      - Lenght of the common suffix in the size
%                                calibration file names.
%   'mostRecentGazeCal'        - Option to identify whether the most recent
%                                GazeCalibration to apply for a given run
%                                is the one 'before' or 'after' the
%                                acquisition of said run (based on files
%                                timestamps).
% Output:
%   None. The routine saves files but does not return variables.
%% Parse input and define variables
p = inputParser; p.KeepUnmatched = true;

% required input
p.addRequired('pathParams',@isstruct);

% optional input
p.addParameter('videoTypeChoice', 'LiveTrackWithVTOP', @ischar);
p.addParameter('variableNamingConvention', 'LiveTrackWithVTOP', @ischar);
p.addParameter('customFunCalls', {}, @iscell);
p.addParameter('skipStageByName', {}, @iscell);
p.addParameter('lastStage', '', @ischar);
p.addParameter('sizeCalIdentifier', '*Scale*_pupil.mat', @ischar);
p.addParameter('sizeCalSuffixLength', 10, @isnumeric);
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
        
        pupilFileNameCAL = fullfile(pathParams.dataOutputDirFull, [GazeCalNameRoot '_pupil.mat']);
        glintFileNameCAL = fullfile(pathParams.dataOutputDirFull, [GazeCalNameRoot '_glint.mat']);
        targetsFileNameCAL = fullfile(pathParams.dataOutputDirFull, [GazeCalNameRoot '_targets.mat']);
        gazeDataFileName = fullfile(pathParams.dataOutputDirFull, [pathParams.runName '_gazeCalData.mat']);
        gazeCalFactorsFileName = fullfile(pathParams.dataOutputDirFull, [pathParams.runName '_gazeCalFactors.mat']);
        
        % for gaze calibration
        calibratedGazeFileName = fullfile(pathParams.dataOutputDirFull, [pathParams.runName '_calibratedGaze.mat']);
        
        % for timebase
        ltReportFileName = fullfile(pathParams.dataSourceDirFull, [pathParams.runName '_report.mat']);
        timebaseFileName = fullfile(pathParams.dataOutputDirFull, [pathParams.runName '_timebase.mat']);
        
    case 'LiveTrackOnlyForGaze'
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
        rawDataPath =  pathParams.dataOutputDirFull;   
        gazeDataFileName = fullfile(pathParams.dataOutputDirFull, [pathParams.runName '_gazeCalData.mat']);
        gazeCalFactorsFileName = fullfile(pathParams.dataOutputDirFull, [pathParams.runName '_gazeCalFactors.mat']);
        
        % for gaze calibration
        calibratedGazeFileName = fullfile(pathParams.dataOutputDirFull, [pathParams.runName '_calibratedGaze.mat']);
        
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
            'applySizeCalibration(pupilFileName,sizeCalFactorsFileName,calibratedPupilFileName, varargin{:});'...
            'makeTargetsFile(LTdatFileName,gazeDataFileName,varargin{:});'...
            'calcGazeCalData(pupilFileNameCAL,glintFileNameCAL,targetsFileNameCAL,gazeDataFileName,varargin{:})' ...
            'calcGazeCalFactors(gazeDataFileName,gazeCalFactorsFileName, varargin{:});' ...
            'applyGazeCalibration(pupilFileName,glintFileName,gazeCalFactorsFileName,calibratedGazeFileName, varargin{:});' ...
            };
    case 'LiveTrackOnlyForGaze'
        funCalls = {...
            'deriveTimebaseFromLTData(glintFileName,ltReportFileName,timebaseFileName, varargin{:});'...
            'calcSizeCalFactors(sizeDataFilesNames, sizeCalFactorsFileName, varargin{:});'...
            'applySizeCalibration(pupilFileName,sizeCalFactorsFileName,calibratedPupilFileName, varargin{:});'...
            'copyLTdatToGazeData(LTdatFileName,gazeDataFileName,varargin{:});'...
            'calcGazeCalFactors(gazeDataFileName,gazeCalFactorsFileName, varargin{:});' ...
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
    if ~any(strcmp(p.Results.skipStageByName,funNames{ff}))
        if strcmp(funNames{ff},'prepareLTGazeCalibrationData') && isempty(LTdatFileName)
            break 
        end
        eval(funCalls{ff});
        if strcmp(p.Results.lastStage,funNames{ff})
            return
        end
    end
end

end % main function


function [LTdatFileName,GazeCalNameRoot] = pickLTGazeData(pathParams,mostRecentGazeCal)
% this function will automatically select the calibration that is most
% appropriate for the run

gazeData = dir(fullfile(pathParams.dataSourceDirFull,'*LTdat*.mat'));
if length(gazeData) == 1
    LTdatFileName = fullfile(pathParams.dataSourceDirFull,gazeData.name);
elseif length(gazeData) > 1
    % sort Gaze Data files by timestamp
    [~,idx] = sort([gazeData.datenum]);
    % get the timestamp for the current run from the corresponding
    % report
    reportFile = dir(fullfile(pathParams.dataSourceDirFull,[pathParams.runName '_report.mat']));
    reportTime = reportFile.datenum;
    switch mostRecentGazeCal
        case 'before'
            % take the most recent calibration file acquired before the current run
            for ii = 1: length(idx)
                if gazeData(idx(ii)).datenum < reportTime
                    LTdatFileName = (fullfile(pathParams.dataSourceDirFull,gazeData(idx(ii)).name));
                    GazeCalNameRoot = gazeData(idx(ii)).name(1:end-10);
                else
                    continue
                end
            end
    end
else
    warning('No gaze calibration found for this run. Skipping gaze calibration step')
    LTdatFileName = '';
end
end

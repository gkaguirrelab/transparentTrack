function processVideoPipeline( pathParams, varargin )
% 
% this is the standard processing pipeline for eye tracking videos. 
% 
% The pipeline consists in the following stages:
%   raw2gray
%   trackGlint
%   extractPupilPerimeter
%   makePreliminaryControlFile
%   correctPupilPerimeter
%   bayesFitPupilPerimeter
%   fitIrisAndPalpebralFissure
%   makePupilFitVideo
% 
% The user can stop the execution after any of the stages with the optional
% param 'lastStage', or skip any amount of stages listing them in a cell
% under the optional param 'skipStage'. Every stage however requires the
% existence of the output from the preceeding ones to be correctly
% executed.

%% Parse input and define variables
p = inputParser; p.KeepUnmatched = true;

% required input
p.addRequired('pathParams',@isstruct);

% optional input
p.addParameter('lastStage', 'makePupilFitVideo', @ischar);
p.addParameter('skipStage', {}, @iscell);

% parse
p.parse(pathParams, varargin{:})
pathParams=p.Results.pathParams;


%% Create output directories if needed
if ~exist(pathParams.dataOutputDirFull,'dir')
    mkdir(pathParams.dataOutputDirFull)
end
if ~exist(pathParams.controlFileDirFull,'dir')
    mkdir(pathParams.controlFileDirFull)
end


%% Define input and output filenames

% Determine if the suffix of the raw file is "_raw.mov" or ".mov"
if exist(fullfile(pathParams.dataSourceDirFull,[pathParams.runName '_raw.mov']),'file')
    rawVideoName = fullfile(pathParams.dataSourceDirFull,[pathParams.runName '_raw.mov']);
else
    rawVideoName = fullfile(pathParams.dataSourceDirFull,[pathParams.runName '.mov']);
end

grayVideoName = fullfile(pathParams.dataOutputDirFull, [pathParams.runName '_gray.avi']);
glintFileName = fullfile(pathParams.dataOutputDirFull, [pathParams.runName '_glint.mat']);
perimeterFileName = fullfile(pathParams.dataOutputDirFull, [pathParams.runName '_perimeter.mat']);
controlFileName = fullfile(pathParams.controlFileDirFull, [pathParams.runName '_controlFile.csv']);
correctedPerimeterFileName = fullfile(pathParams.dataOutputDirFull, [pathParams.runName '_correctedPerimeter.mat']);
ellipseFitFileName = fullfile(pathParams.dataOutputDirFull, [pathParams.runName '_pupil.mat']);
finalFitVideoName = fullfile(pathParams.dataOutputDirFull, [pathParams.runName '_finalFit.mat']);
irisFitFileName = fullfile(pathParams.dataOutputDirFull, [pathParams.runName '_iris.mat']);
palpebralFissureFileName = fullfile(pathParams.dataOutputDirFull, [pathParams.runName '_palpebralFissure.mat']);

%% Conduct the analysis
% NOTE: some of the analysis steps are wrapped in a while+try/catch loop to
% fix the matlabprefs.mat corruption bug.

% Convert raw video to cropped, resized, 60Hz gray
if ~any(strcmp(p.Results.skipStage,'raw2gray'))
    raw2gray(rawVideoName,grayVideoName, varargin{:});
    if strcmp(p.Results.lastStage,'raw2gray')
        return
    end
end

% track the glint
if ~any(strcmp(p.Results.skipStage,'trackGlint'))
    % intialize while control
    success = 0;
    while ~success
        try
            trackGlint(grayVideoName, glintFileName, varargin{:});
            success = 1;
        catch ME
            % if there is a corruption error clear matlabprefs.mat and try again
            if regexp(ME.message, ...
                    regexptranslate('*matlabprefs.mat. File might be corrupt.'))
                warning ('File matlabprefs.mat corrupt during execution. Cleaning up and trying again.')
                matlabprefsCleanup;
                success = 0;
            elseif regexp(ME.message, ...
                    regexptranslate('Error detected on workers*'))
                warning ('Error detected on workers. Cleaning up matlabprefs and trying again.')
                matlabprefsCleanup;
                success = 0;
            % if a parpool is already open, close it and try again
            elseif (strcmp(ME.message, ...
                    'Found an interactive session. You cannot have multiple interactive sessions open simultaneously. To terminate the existing session, use ''delete(gcp(''nocreate''))''.'))
                warning ('Found a parpool already open. Closing it and trying again.')
                delete(gcp('nocreate'))
                success = 0;
            else
                rethrow(ME)
            end
        end
    end
    % clear while control
    clear success
    if strcmp(p.Results.lastStage,'trackGlint')
        return
    end
end

% extract pupil perimeter
if ~any(strcmp(p.Results.skipStage,'extractPupilPerimeter'))
    % intialize while control
    success = 0;
    while ~success
        try
            extractPupilPerimeter(grayVideoName, perimeterFileName, varargin{:});
            success = 1;
        catch ME
            % if there is a corruption error clear matlabprefs.mat and try again
            if regexp(ME.message, ...
                    regexptranslate('*matlabprefs.mat. File might be corrupt.'))
                warning ('File matlabprefs.mat corrupt during execution. Cleaning up and trying again.')
                matlabprefsCleanup;
                success = 0;
            elseif regexp(ME.message, ...
                    regexptranslate('Error detected on workers*'))
                warning ('Error detected on workers. Cleaning up matlabprefs and trying again.')
                matlabprefsCleanup;
                success = 0;
            % if a parpool is already open, close it and try again
            elseif (strcmp(ME.message, ...
                    'Found an interactive session. You cannot have multiple interactive sessions open simultaneously. To terminate the existing session, use ''delete(gcp(''nocreate''))''.'))
                warning ('Found a parpool already open. Closing it and trying again.')
                delete(gcp('nocreate'))
                success = 0;
            else
                rethrow(ME)
            end
        end
    end
    % clear while control
    clear success
    if strcmp(p.Results.lastStage,'extractPupilPerimeter')
        return
    end
end

% generate preliminary control file
if ~any(strcmp(p.Results.skipStage,'makePreliminaryControlFile'))
    % intialize while control
    success = 0;
    while ~success
        try
            makePreliminaryControlFile(controlFileName, perimeterFileName, glintFileName, varargin{:});
            success = 1;
        catch ME
            % if there is a corruption error clear matlabprefs.mat and try again
            if regexp(ME.message, ...
                    regexptranslate('*matlabprefs.mat. File might be corrupt.'))
                warning ('File matlabprefs.mat corrupt during execution. Cleaning up and trying again.')
                matlabprefsCleanup;
                success = 0;
            elseif regexp(ME.message, ...
                    regexptranslate('Error detected on workers*'))
                warning ('Error detected on workers. Cleaning up matlabprefs and trying again.')
                matlabprefsCleanup;
                success = 0;
            % if a parpool is already open, close it and try again
            elseif (strcmp(ME.message, ...
                    'Found an interactive session. You cannot have multiple interactive sessions open simultaneously. To terminate the existing session, use ''delete(gcp(''nocreate''))''.'))
                warning ('Found a parpool already open. Closing it and trying again.')
                delete(gcp('nocreate'))
                success = 0;
            else
                rethrow(ME)
            end
        end
    end
    % clear while control
    clear success
    if strcmp(p.Results.lastStage,'makePreliminaryControlFile')
        return
    end
end

% correct the perimeter video
if ~any(strcmp(p.Results.skipStage,'correctPupilPerimeter'))
    % intialize while control
    success = 0;
    while ~success
        try
            correctPupilPerimeter(perimeterFileName,controlFileName,correctedPerimeterFileName, varargin{:});
            success = 1;
        catch ME
            % if there is a corruption error clear matlabprefs.mat and try again
            if regexp(ME.message, ...
                    regexptranslate('*matlabprefs.mat. File might be corrupt.'))
                warning ('File matlabprefs.mat corrupt during execution. Cleaning up and trying again.')
                matlabprefsCleanup;
                success = 0;
            elseif regexp(ME.message, ...
                    regexptranslate('Error detected on workers*'))
                warning ('Error detected on workers. Cleaning up matlabprefs and trying again.')
                matlabprefsCleanup;
                success = 0;
            % if a parpool is already open, close it and try again
            elseif (strcmp(ME.message, ...
                    'Found an interactive session. You cannot have multiple interactive sessions open simultaneously. To terminate the existing session, use ''delete(gcp(''nocreate''))''.'))
                warning ('Found a parpool already open. Closing it and trying again.')
                delete(gcp('nocreate'))
                success = 0;
            else
                rethrow(ME)
            end
        end
    end
    % clear while control
    clear success
    if strcmp(p.Results.lastStage,'correctPupilPerimeter')
        return
    end
end

% bayesian fit of the pupil on the corrected perimeter video
if ~any(strcmp(p.Results.skipStage,'bayesFitPupilPerimeter'))
    % intialize while control
    success = 0;
    while ~success
        try
            bayesFitPupilPerimeter(correctedPerimeterFileName, ellipseFitFileName, varargin{:});
            success = 1;
        catch ME
            % if there is a corruption error clear matlabprefs.mat and try again
            if regexp(ME.message, ...
                    regexptranslate('*matlabprefs.mat. File might be corrupt.'))
                warning ('File matlabprefs.mat corrupt during execution. Cleaning up and trying again.')
                matlabprefsCleanup;
                success = 0;
            elseif regexp(ME.message, ...
                    regexptranslate('Error detected on workers*'))
                warning ('Error detected on workers. Cleaning up matlabprefs and trying again.')
                matlabprefsCleanup;
                success = 0;
            % if a parpool is already open, close it and try again
            elseif (strcmp(ME.message, ...
                    'Found an interactive session. You cannot have multiple interactive sessions open simultaneously. To terminate the existing session, use ''delete(gcp(''nocreate''))''.'))
                warning ('Found a parpool already open. Closing it and trying again.')
                delete(gcp('nocreate'))
                success = 0;
            else
                rethrow(ME)
            end
        end
    end
    % clear while control
    clear success
    if strcmp(p.Results.lastStage,'bayesFitPupilPerimeter')
        return
    end
end

% fit Iris and palpebral fissure
if ~any(strcmp(p.Results.skipStage,'fitIrisAndPalpebralFissure'))
    % intialize while control
    success = 0;
    while ~success
        try
            fitIrisAndPalpebralFissure(grayVideoName, perimeterFileName, ellipseFitFileName, irisFitFileName, palpebralFissureFileName, varargin{:});
            success = 1;
        catch ME
            % if there is a corruption error clear matlabprefs.mat and try again
            if regexp(ME.message, ...
                    regexptranslate('*matlabprefs.mat. File might be corrupt.'))
                warning ('File matlabprefs.mat corrupt during execution. Cleaning up and trying again.')
                matlabprefsCleanup;
                success = 0;
            elseif regexp(ME.message, ...
                    regexptranslate('Error detected on workers*'))
                warning ('Error detected on workers. Cleaning up matlabprefs and trying again.')
                matlabprefsCleanup;
                success = 0;
            % if a parpool is already open, close it and try again
            elseif (strcmp(ME.message, ...
                    'Found an interactive session. You cannot have multiple interactive sessions open simultaneously. To terminate the existing session, use ''delete(gcp(''nocreate''))''.'))
                warning ('Found a parpool already open. Closing it and trying again.')
                delete(gcp('nocreate'))
                success = 0;
            else
                rethrow(ME)
            end
        end
    end
    if strcmp(p.Results.lastStage,'fitIrisAndPalpebralFissure')
        return
    end
end

% create a video of the final fit
if ~any(strcmp(p.Results.skipStage,'makePupilFitVideo'))
    % intialize while control
    success = 0;
    while ~success
        try
            makePupilFitVideo(grayVideoName, finalFitVideoName, ...
        'glintFileName', glintFileName, 'perimeterFileName', correctedPerimeterFileName,...
        'ellipseFitFileName', ellipseFitFileName, 'whichFieldToPlot', 'pPosteriorMeanTransparent', ...
        'irisFitFileName', irisFitFileName, ...
        'controlFileName',controlFileName,varargin{:});
            success = 1;
        catch ME
            % if there is a corruption error clear matlabprefs.mat and try again
            if regexp(ME.message, ...
                    regexptranslate('*matlabprefs.mat. File might be corrupt.'))
                warning ('File matlabprefs.mat corrupt during execution. Cleaning up and trying again.')
                matlabprefsCleanup;
                success = 0;
            elseif regexp(ME.message, ...
                    regexptranslate('Error detected on workers*'))
                warning ('Error detected on workers. Cleaning up matlabprefs and trying again.')
                matlabprefsCleanup;
                success = 0;
            % if a parpool is already open, close it and try again
            elseif (strcmp(ME.message, ...
                    'Found an interactive session. You cannot have multiple interactive sessions open simultaneously. To terminate the existing session, use ''delete(gcp(''nocreate''))''.'))
                warning ('Found a parpool already open. Closing it and trying again.')
                delete(gcp('nocreate'))
                success = 0;
            else
                rethrow(ME)
            end
        end
    end
    % clear while control
    clear success
end

end % function


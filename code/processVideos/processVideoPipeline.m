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
pupilFileName = fullfile(pathParams.dataOutputDirFull, [pathParams.runName '_pupil.mat']);
finalFitVideoName = fullfile(pathParams.dataOutputDirFull, [pathParams.runName '_finalFit.avi']);
irisFileName = fullfile(pathParams.dataOutputDirFull, [pathParams.runName '_iris.mat']);

%% Conduct the analysis
% NOTE: some of the analysis steps are wrapped in a while+try/catch loop to
% fix the matlabprefs.mat corruption bug.

funCalls = {...
    'raw2gray(rawVideoName,grayVideoName, varargin{:});'...
    'trackGlint(grayVideoName, glintFileName, varargin{:});'...
    'trackGlint(grayVideoName, glintFileName, varargin{:});'...
    'extractPupilPerimeter(grayVideoName, perimeterFileName, varargin{:});'...
    'correctPupilPerimeter(perimeterFileName,controlFileName,correctedPerimeterFileName, varargin{:});' ...
    'bayesFitPupilPerimeter(correctedPerimeterFileName, pupilFileName, varargin{:});'...
    'fitIrisCircleAndMask(grayVideoName, perimeterFileName, pupilFileName, irisFileName, varargin{:});' ...
    ['makePupilFitVideo(grayVideoName, finalFitVideoName,' ...
    '''glintFileName'', glintFileName, ''perimeterFileName'', correctedPerimeterFileName,'...
    '''pupilFileName'', pupilFileName, ''whichFieldToPlot'', ''pPosteriorMeanTransparent'',' ...
    '''irisFileName'', irisFileName,' ...
    '''controlFileName'',controlFileName,varargin{:});']...
    };

funNames = { ...
    'raw2gray'...
    'trackGlint'...
    'extractPupilPerimeter'...
    'makePreliminaryControlFile'...
    'correctPupilPerimeter'...
    'bayesFitPupilPerimeter'...
    'fitIrisAndPalpebralFissure'...
    'makePupilFitVideo'...
    };

for ff = 1:length(funCalls)
    if ~any(strcmp(p.Results.skipStage,funNames{ff}))
        % intialize while control
        success = 0;
        attempts = 1; % we will attempt to execute the instruction 3 times. If it does not work, the code will eventually break.
        while ~success
            try
                eval(funCalls{ff});
                success = 1;
            catch ME
                warning ('There has been an error during execution. Cleaning matlabprefs.mat and trying again - attempt %d.',attempts)
                matlabprefsCleanup;
                success = 0;
                attempts = attempts +1;
                if attempts > 3
                    rethrow(ME)
                end
            end
        end
        % clear while and if control
        clear success
        clear attempts
        if strcmp(p.Results.lastStage,funNames{ff})
            return
        end
    end
end
end

function processVideoPipeline( pathParams, varargin )

%% Parse input and define variables
p = inputParser; p.KeepUnmatched = true;

% required input
p.addRequired('pathParams',@isstruct);

% optional input
p.addParameter('lastStage', 'makePupilFitVideo', @ischar);


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


%% Define inout and output filenames

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


%% Conduct the analysis

% Convert raw video to cropped, resized, 60Hz gray
raw2gray(rawVideoName,grayVideoName, varargin{:});
if strcmp(p.Results.lastStage,'raw2gray')
    return
end

% track the glint
trackGlint(grayVideoName, glintFileName, varargin{:});
if strcmp(p.Results.lastStage,'trackGlint')
    return
end

% extract pupil perimeter
extractPupilPerimeter(grayVideoName, perimeterFileName, varargin{:});
if strcmp(p.Results.lastStage,'extractPupilPerimeter')
    return
end

% generate preliminary control file
makePreliminaryControlFile(controlFileName, perimeterFileName, glintFileName, varargin{:});
if strcmp(p.Results.lastStage,'makePreliminaryControlFile')
    return
end

% correct the perimeter video
correctPupilPerimeter(perimeterFileName,controlFileName,correctedPerimeterFileName, varargin{:});
if strcmp(p.Results.lastStage,'correctPupilPerimeter')
    return
end

% bayesian fit of the pupil on the corrected perimeter video
bayesFitPupilPerimeter(correctedPerimeterFileName, ellipseFitFileName, varargin{:});
if strcmp(p.Results.lastStage,'bayesFitPupilPerimeter')
    return
end

% create a video of the final fit
makePupilFitVideo(grayVideoName, finalFitVideoName, ...
    'glintFileName', glintFileName, 'perimeterFileName', correctedPerimeterFileName,...
    'ellipseFitFileName', ellipseFitFileName, 'whichFieldToPlot', 'pPosteriorMeanTransparent', ...
    'controlFileName',controlFileName,varargin{:});

end % function


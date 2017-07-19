function processVideoPipeline( pathParams, varargin )

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
irisFitFileName = fullfile(pathParams.dataOutputDirFull, [pathParams.runName '_iris.mat']);
palpebralFissureFileName = fullfile(pathParams.dataOutputDirFull, [pathParams.runName '_palpebralFissure.mat']);
finalFitVideoName = fullfile(pathParams.dataOutputDirFull, [pathParams.runName '_finalFit.mat']);


%% Conduct the analysis
% 
% % Convert raw video to cropped, resized, 60Hz gray
% if ~any(strcmp(p.Results.skipStage,'raw2gray'))
%     raw2gray(rawVideoName,grayVideoName, varargin{:});
%     if strcmp(p.Results.lastStage,'raw2gray')
%         return
%     end
% end
% 
% % track the glint
% if ~any(strcmp(p.Results.skipStage,'trackGlint'))
%     trackGlint(grayVideoName, glintFileName, varargin{:});
%     if strcmp(p.Results.lastStage,'trackGlint')
%         return
%     end
% end
% 
% % extract pupil perimeter
% if ~any(strcmp(p.Results.skipStage,'extractPupilPerimeter'))
%     extractPupilPerimeter(grayVideoName, perimeterFileName, varargin{:});
%     if strcmp(p.Results.lastStage,'extractPupilPerimeter')
%         return
%     end
% end
% 
% % generate preliminary control file
% if ~any(strcmp(p.Results.skipStage,'makePreliminaryControlFile'))
%     makePreliminaryControlFile(controlFileName, perimeterFileName, glintFileName, varargin{:});
%     if strcmp(p.Results.lastStage,'makePreliminaryControlFile')
%         return
%     end
% end
% 
% % correct the perimeter video
% if ~any(strcmp(p.Results.skipStage,'correctPupilPerimeter'))
%     correctPupilPerimeter(perimeterFileName,controlFileName,correctedPerimeterFileName, varargin{:});
%     if strcmp(p.Results.lastStage,'correctPupilPerimeter')
%         return
%     end
% end
% 
% % bayesian fit of the pupil on the corrected perimeter video
% if ~any(strcmp(p.Results.skipStage,'bayesFitPupilPerimeter'))
%     bayesFitPupilPerimeter(correctedPerimeterFileName, ellipseFitFileName, varargin{:});
%     if strcmp(p.Results.lastStage,'bayesFitPupilPerimeter')
%         return
%     end
% end

% % fit Iris and palpebral fissure
% if ~any(strcmp(p.Results.skipStage,'fitIrisAndPalpebralFissure'))
%     fitIrisAndPalpebralFissure(grayVideoName, perimeterFileName, ellipseFitFileName, irisFitFileName, palpebralFissureFileName, varargin{:});
%     if strcmp(p.Results.lastStage,'fitIrisAndPalpebralFissure')
%         return
%     end
% end

% create a video of the final fit
if ~any(strcmp(p.Results.skipStage,'makePupilFitVideo'))
    makePupilFitVideo(grayVideoName, finalFitVideoName, ...
        'glintFileName', glintFileName, 'perimeterFileName', correctedPerimeterFileName,...
        'ellipseFitFileName', ellipseFitFileName, 'whichFieldToPlot', 'pPosteriorMeanTransparent', ...
        'irisFitFileName', irisFitFileName, ...
        'controlFileName',controlFileName,varargin{:});
end

end % function


function runVideoPipeline( pathParams, varargin )
% runVideoPipeline( pathParams, varargin ) - A standard processing pipeline for eye tracking videos.
%
% The pipeline consists of the following stages:
%   deinterlaceVideo
%   resizeAndCropVideo
%   findGlint
%   findPupilPerimeter
%   makeControlFile
%   applyControlFile
%   fitPupilPerimeter
%   estimateSceneGeometry
%   smoothPupilArea
%   makeFitVideo
%
% The user can stop the execution after any of the stages with the optional
% param 'lastStage', or skip any number of stages listing them in a cell
% under the optional param 'skipStageByName'. Every stage, however, requires the
% existence of the output from the preceeding ones to be correctly
% executed.
%
% INPUT
%   pathParams - This structure has fields corresponding to the name and
%       location of the files to be processed. Fields include:
%           dataSourceDirFull - full path to the directory that contains
%              the source file
%           dataOutputDirFull - full path to the directory where the files
%               that result from processing will be written
%           runName - the stem name for the source and subsequent output
%               files.
%       This structure may have other fields defined that are used in
%       routines prior to this one to assemble the full paths to the source
%       and output dirs.
%
% OUTPUT
%   None. The routine saves files but does not return variables.
%
% OPTIONS
%   lastStage - the last stage to be executed. By deafult ends with the
%      production of the fit video.
%   skipStageByName - a cell array of function calls to be skipped during
%      execution of the pipeline.
%   skipStageByNumber - an array of function calls numbers to be skipped
%      during execution of the pipeline.
%   displayAvailableStages - displays numbered list of the available
%       stages, to be used as reference for skupStageByNumber. Note that if
%       this is set to true, no analysis will be perfomed.
%   rawVideoSuffix - cell array of strings that contain possible suffixes
%       of raw video files to be processed. The routine will search for
%       raw videos sequentially in the cell array until it finds a match
%   videoTypeChoice - This key-value can be used to identify a set of
%       processing choices. There are three choices defined here that are
%       somewhat idiosyncratic to the data being collected in the GKAguirre
%       lab of the University of Pennsylvania. Set this value to 'custom'
%       to execute a set of parameters passed using customFunCalls.
%           LiveTrackWithVTOP_sizeCal - used for the analysis of size
%           calibration videos of a black circle on a calibration wand.
%           Skips tracking of the glint or production of a control file, as
%           no eyelid will be present.
%           LiveTrackWithVTOP_eye - our standard processing pipeline.
%           LiveTrackWithVTOP_eyeNoIris - the standard, but skipping
%           fitting of the iris bounds as this is under development.
%   customFunCalls - A cell array of functional calls and key values that
%       can be passed in lieu of hard-coding a videoTypeChoice set here.
%       This is used in concert with passing 'custom' to videoTypeChoice.
%   catchErrors - controls if the function calls take place within a
%       try-catch block. If set to true (the default) then the routine will
%       attempt to execute a function three times before exiting with an
%       error. After each error within the try-catch block, the function
%       cleanupMatlabPrefs is called. This is thought to correct a
%       stochastic error that can occur in parpool jobs and results in the
%       corruption of the matlab preference file, which is then deleted.
%    maxAttempts - the number of times that a given stage will be re-tried
%       in the event of an error.
%    makeFitVideoByName - a cell array of stages for which a fit
%       video will be produced following completion of the stage. The video
%       file name is the run name, followed by "_fitStageX.avi" where X is
%       the idx of the function call list.
%    makeFitVideoByNumber - a cell array of stages for which a fit
%       video will be produced following completion of the stage. The video
%       file name is the run name, followed by "_fitStageX.avi" where X is
%       the idx of the function call list.
%


%% Parse input and define variables
p = inputParser; p.KeepUnmatched = true;

% required input
p.addRequired('pathParams',@isstruct);

% optional input
p.addParameter('lastStage', 'makePupilFitVideo', @ischar);
p.addParameter('skipStageByName', {}, @iscell);
p.addParameter('skipStageByNumber',[], @isnumeric);
p.addParameter('displayAvailableStages', false, @islogical)
p.addParameter('rawVideoSuffix', {'_raw.mov' '.mov'}, @iscell);
p.addParameter('videoTypeChoice', 'LiveTrackWithVTOP_eye', @ischar);
p.addParameter('customFunCalls', {}, @iscell);
p.addParameter('catchErrors', true, @islogical);
p.addParameter('maxAttempts',3,@isnueric);
p.addParameter('makeFitVideoByName',{},@iscell);
p.addParameter('makeFitVideoByNumber',[],@isnumeric);

% parse
p.parse(pathParams, varargin{:})
pathParams=p.Results.pathParams;

%% Sanity check the input
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


%% Define input filenames
if ~any(strcmp(p.Results.skipStageByName,'deinterlaceVideo')) && ~any(p.Results.skipStageByNumber == 1) 
    % Create a cell array of candidate raw video nmaes with the runName and
    % each of the rawVideoSuffix choices
    candidateRawVideoNames = ...
        cellfun(@(x) fullfile(pathParams.dataSourceDirFull,[pathParams.runName x]),p.Results.rawVideoSuffix,'uniformoutput',false);
    
    % Test if each of these candidate video files exist
    existTest = cellfun(@(x) exist(x,'file')==2, candidateRawVideoNames);
    
    switch (sum(existTest))
        case 0
            error('Cannot find a raw video with this run name and the specified suffix');
        case 1
            rawVideoName = candidateRawVideoNames{ existTest==1 };
        otherwise
            error('There is more than one raw video with this run name and the specified suffix');
    end
end

%% Define output filenames
grayVideoName = fullfile(pathParams.dataOutputDirFull, [pathParams.runName '_gray.avi']);
glintFileName = fullfile(pathParams.dataOutputDirFull, [pathParams.runName '_glint.mat']);
perimeterFileName = fullfile(pathParams.dataOutputDirFull, [pathParams.runName '_perimeter.mat']);
controlFileName = fullfile(pathParams.dataOutputDirFull, [pathParams.runName '_controlFile.csv']);
correctedPerimeterFileName = fullfile(pathParams.dataOutputDirFull, [pathParams.runName '_correctedPerimeter.mat']);
pupilFileName = fullfile(pathParams.dataOutputDirFull, [pathParams.runName '_pupil.mat']);
sceneGeometryFileName = fullfile(pathParams.dataOutputDirFull, [pathParams.runName '_sceneGeometry.mat']);
sceneDiagnosticPlotFileName = fullfile(pathParams.dataOutputDirFull, [pathParams.runName '_sceneDiagnosticPlot.pdf']);
finalFitVideoName = fullfile(pathParams.dataOutputDirFull, [pathParams.runName '_finalFit.avi']);
irisFileName = fullfile(pathParams.dataOutputDirFull, [pathParams.runName '_iris.mat']);


%% Conduct the analysis
% We define a set of function calls, and then execute them within a try-
%  catch block. This is chiefly because a Matlab bug causes the matlab
%  prefs file to sometimes become corrupted during parpool operations.
%  If we detect an error, we delete the pref file and try again.

% define analysis pipelines.
switch p.Results.videoTypeChoice
    case 'LiveTrackWithVTOP_sizeCal'
        funCalls = {...
            'deinterlaceVideo(rawVideoName,grayVideoName, varargin{:});'...
            'findPupilPerimeter(grayVideoName, perimeterFileName, varargin{:});'...
            ['fitPupilPerimeter(perimeterFileName, pupilFileName,' ...
            '''nSplits'', 0, varargin{:});']...
            ['makeFitVideo(grayVideoName, finalFitVideoName,' ...
            ' ''perimeterFileName'', perimeterFileName,''perimeterColor'',''r'','...
            ' ''pupilFileName'', pupilFileName, ''whichFieldToPlot'', ''ellipseParamsUnconstrained_mean'',' ...
            ' varargin{:});']...
            };
    case 'LiveTrackWithVTOP_eye'
        funCalls = {...
            'deinterlaceVideo(rawVideoName,grayVideoName, varargin{:});'...
            'findGlint(grayVideoName, glintFileName, varargin{:});'...
            'findPupilPerimeter(grayVideoName, perimeterFileName, varargin{:});'...
            'makeControlFile(controlFileName, perimeterFileName, glintFileName, varargin{:});' ...
            'applyControlFile(perimeterFileName,controlFileName,correctedPerimeterFileName, varargin{:});' ...
            ['fitPupilPerimeter(correctedPerimeterFileName, pupilFileName,' ...
            '''nSplits'', 0, varargin{:});']...
            ['estimateSceneGeometry(pupilFileName, sceneGeometryFileName,' ...
            '''sceneDiagnosticPlotFileName'', sceneDiagnosticPlotFileName, varargin{:});']...
            ['fitPupilPerimeter(correctedPerimeterFileName, pupilFileName,' ...
            '''sceneGeometryFileName'', sceneGeometryFileName, varargin{:});']...
            'smoothPupilArea(correctedPerimeterFileName, pupilFileName, sceneGeometryFileName, varargin{:});'...
            ['makeFitVideo(grayVideoName, finalFitVideoName,' ...
            '''glintFileName'', glintFileName, ''perimeterFileName'', correctedPerimeterFileName,'...
            '''pupilFileName'', pupilFileName, ''sceneGeometryFileName'', sceneGeometryFileName,' ...
            '''controlFileName'',controlFileName,varargin{:});']...
            };
    case 'custom'
        funCalls = p.Results.customFunCalls;
    otherwise
        error('VideoTypeChoice is not defined')
end

% Grab the function names as the portion of the funCalls that preceed the
% open parenthesis
funNames = cellfun(@(x) strtok(x,'('),funCalls,'UniformOutput',false);

% check if just displaying stages in use
if p.Results.displayAvailableStages
    display(['Available stages for ' p.Results.videoTypeChoice])
    for ff = 1:length(funCalls)
        fprintf ([num2str(ff) ' - ' funNames{ff} '\n'])
    end
    return
end

% Loop through the function calls
for ff = 1:length(funCalls)
    
    % Check if we are instructed to skip this stage, either by stage name
    % or by functional call number.
    if ~any(strcmp(p.Results.skipStageByName,funNames{ff})) && ~any(p.Results.skipStageByNumber == ff)
        % check if we are going to catch or throw errors
        if p.Results.catchErrors
            success = false; % flag for successful stage execution
            attempts = 1; % which attempt are we on?
            % intialize while control
            while ~success
                try
                    eval(funCalls{ff});
                    success = true;
                catch ME
                    warning ('There has been an error during execution. Cleaning matlabprefs.mat and trying again - attempt %d.',attempts)
                    cleanupMatlabPrefs;
                    success = false;
                    attempts = attempts +1;
                    if attempts > p.Results.maxAttempts
                        rethrow(ME)
                    end
                end
            end
            % clear while and if control
            clear success
            clear attempts
        else % we will throw errors
            eval(funCalls{ff});
        end % if catchErrors
        
        % Check if we should make a fit video for this stage
        if any(strcmp(p.Results.makeFitVideoByName,funNames{ff})) || any(p.Results.makeFitVideoByNumber == ff)
            makeFitVideoForThisStage(pathParams, funNames, ff, varargin{:});
        end
        
        % clear all files (hopefully prevents 'too many files open' error)
        fclose all;
        if strcmp(p.Results.lastStage,funNames{ff})
            break
        end
    end % if we aren't skipping this stage by name or number
end % loop over function calls


end % main function


%% LOCAL FUNCTIONS

function cleanupMatlabPrefs
% This small function closes any open parallel pool and removes the
% matlabprefs.mat file in case it is corrupt.
% This file can get corrupt while using parallel pools because all of the
% pool worker try to read/write it at the same time. MathWorks does not
% currently have a better solution for this probem:
%   https://www.mathworks.com/matlabcentral/answers/348966-how-to-set-prefdir-for-parpool-workers

% close parallel pool if any is open
poolobj = gcp('nocreate');
if ~isempty (poolobj)
    delete(poolobj);
end

% clean up the corrupt matlabprefs file
prefFile = fullfile(prefdir,'matlabprefs.mat');
system(['rm -rf "' prefFile '"'])
end % cleanupMatlabPrefs


function makeFitVideoForThisStage(pathParams, funNames, ff, varargin)

% Define the fitVideo output name
fitVideoFileName = fullfile(pathParams.dataOutputDirFull, [pathParams.runName '_fitStage' num2str(ff) '.avi']);

% Assemble the entire list of potential files to include in the video
grayVideoName = fullfile(pathParams.dataOutputDirFull, [pathParams.runName '_gray.avi']);
glintFileName = fullfile(pathParams.dataOutputDirFull, [pathParams.runName '_glint.mat']);
initialPerimeterFileName = fullfile(pathParams.dataOutputDirFull, [pathParams.runName '_perimeter.mat']);
controlFileName = fullfile(pathParams.dataOutputDirFull, [pathParams.runName '_controlFile.csv']);
correctedPerimeterFileName = fullfile(pathParams.dataOutputDirFull, [pathParams.runName '_correctedPerimeter.mat']);
pupilFileName = fullfile(pathParams.dataOutputDirFull, [pathParams.runName '_pupil.mat']);
sceneGeometryFileName = fullfile(pathParams.dataOutputDirFull, [pathParams.runName '_sceneGeometry.mat']);
irisFileName = fullfile(pathParams.dataOutputDirFull, [pathParams.runName '_iris.mat']);

% Depending upon which stage just completed, set to empty any files
% that do not yet exist, and select the appropriate perimeter file and
% pupil ellipse fit to display.
switch funNames{ff}
    case 'deinterlaceVideo'
        glintFileName = []; perimeterFileName=[]; controlFileName=[]; pupilFileName=[]; sceneGeometryFileName=[]; irisFileName=[];
    case 'resizeAndCropVideo'
        glintFileName = []; perimeterFileName=[]; controlFileName=[]; pupilFileName=[]; sceneGeometryFileName=[]; irisFileName=[];
    case 'findGlint'
        perimeterFileName=[]; controlFileName=[]; pupilFileName=[]; sceneGeometryFileName=[]; irisFileName=[];
    case 'findPupilPerimeter'
        if ~exist(glintFileName,'file')
            glintFileName = [];
        end
        perimeterFileName=initialPerimeterFileName;
        controlFileName=[]; pupilFileName=[]; sceneGeometryFileName=[]; irisFileName=[];
    case 'makeControlFile'
        if ~exist(glintFileName,'file')
            glintFileName = [];
        end
        perimeterFileName=initialPerimeterFileName;
        controlFileName=[]; pupilFileName=[]; sceneGeometryFileName=[]; irisFileName=[];
    case 'applyControlFile'
        if ~exist(glintFileName,'file')
            glintFileName = [];
        end
        perimeterFileName=correctedPerimeterFileName;
        pupilFileName=[]; sceneGeometryFileName=[]; irisFileName=[];
    case 'fitPupilPerimeter'
        if ~exist(glintFileName,'file')
            glintFileName = [];
        end
        perimeterFileName=correctedPerimeterFileName;
        sceneGeometryFileName=[]; irisFileName=[];
        % If the sceneGeometry has been determined by this point, we can
        % plot the sceneConstrained ellipse fit, otherwise plot the
        % unconstrained
        varargin={varargin{:}, 'whichFieldToPlot', 'ellipseParamsUnconstrained_mean'};
        sceneFunCallIdx=find(strcmp(funNames,'estimateSceneGeometry'));
        if ~isempty(sceneFunCallIdx)
            if sceneFunCallIdx < ff
                varargin={varargin{:}, 'whichFieldToPlot', 'ellipseParamsSceneConstrained_mean'};
            end
        end
    case 'estimateSceneGeometry'
        if ~exist(glintFileName,'file')
            glintFileName = [];
        end
        perimeterFileName=correctedPerimeterFileName;
        irisFileName=[];
        varargin={varargin{:}, 'whichFieldToPlot', 'ellipseParamsUnconstrained_mean'};
    case 'smoothPupilArea'
        if ~exist(glintFileName,'file')
            glintFileName = [];
        end
        perimeterFileName=correctedPerimeterFileName;
        irisFileName=[];
    case 'fitIrisPerimeter'
        if ~exist(glintFileName,'file')
            glintFileName = [];
        end
        perimeterFileName=correctedPerimeterFileName;
    otherwise
        warning('I do not recognize that stage. Returning without creating a stage fit video');
        return
end

% make the video
makeFitVideo(grayVideoName, fitVideoFileName, ...
    'glintFileName', glintFileName, 'perimeterFileName', perimeterFileName,...
    'controlFileName',controlFileName, 'pupilFileName', pupilFileName, ...
    'sceneGeometryFileName', sceneGeometryFileName, 'irisFileName', irisFileName, varargin{:});

end % makeFitVideoForThisStage




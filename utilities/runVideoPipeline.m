function runVideoPipeline( pathParams, varargin )
% A standard processing pipeline for eye tracking videos.
%
% Syntax:
%  runVideoPipeline( pathParams )
%
% Description:
%	The pipeline consists of the following stages:
%       deinterlaceVideo
%       findGlint
%       findPupilPerimeter
%       makeControlFile
%       applyControlFile
%       fitPupilPerimeter -- with minimal constraints
%       estimateSceneParams
%       fitPupilPerimeter -- fit again with scene geometry constraints
%       smoothPupilRadius
%       makeFitVideo
%       makeEyeModelVideo
%
%   The user can stop the execution after any of the stages with the
%   optional param 'lastStage', or skip any number of stages listing them
%   in a cell under the optional param 'skipStageByName'. Every stage,
%   however, requires the existence of the output from the preceeding ones
%   to be correctly executed.
%
% Input:
%	pathParams            - This structure has fields corresponding to the
%                           name and location of the files to be processed.
%                           Fields include:
%                           dataSourceDirFull: full path to the directory
%                               that contains the source file
%                           dataOutputDirFull: full path to the directory
%                               where the files that result from processing
%                               will be written
%                           runName: the stem name for the source and
%                               subsequent output files.
%                           This structure may have other fields defined
%                           that are used in routines prior to this one to
%                           assemble the full paths to the source and
%                           output dirs.
%
% Optional key/value pairs (display and I/O):
%  'lastStageByName'      - The last stage to be executed. By deafult ends 
%                           with the production of the fit video.
%  'skipStageByName'      - A cell array of function calls to be skipped 
%                           during execution of the pipeline.
%  'skipStageByNumber'    - An array of function calls numbers to be 
%                           skipped during execution of the pipeline.
%  'displayAvailableStages' - Displays numbered list of the available
%                           stages, to be used as reference for
%                           skipStageByNumber. Note that if this is set to
%                           true, no analysis will be perfomed.
%  'rawVideoSuffix'       - Cell array of strings that contain possible
%                           suffixes of raw video files to be processed.
%                           The routine will search for raw videos
%                           sequentially in the cell array until it finds a
%                           match
%  'videoTypeChoice'      - This key-value can be used to identify a set of
%                           processing choices. There are two choices
%                           defined here that are somewhat idiosyncratic to
%                           the data being collected in the GKAguirre lab
%                           of the University of Pennsylvania. Set this
%                           value to 'custom' to execute a set of
%                           parameters passed using customFunCalls.
%                               LiveTrackWithVTOP_sizeCal: used for the
%                               analysis of size calibration videos of a
%                               black circle on a calibration wand. Skips
%                               tracking of the glint or production of a
%                               control file, as no eyelid will be present
%                               LiveTrackWithVTOP_eye: standard pipeline
%  'customFunCalls'       - A cell array of functional calls and key values
%                           that can be passed in lieu of hard-coding a
%                           videoTypeChoice set here. This is used in
%                           concert with passing 'custom' to
%                           videoTypeChoice.
%  'customSceneGeometryFile' - When passed, this full path to a
%                           sceneGeometry file is passed as input to stages
%                           subsequent to estimateSceneParams. This
%                           allows processing upon one acquisition to use
%                           the scene Geometry derived from a different
%                           run.
%  'catchErrors'          - Controls if the function calls take place 
%                           within a try-catch block. If set to true (the
%                           default) then the routine will attempt to
%                           execute a function 'maxAttempts' times before
%                           exiting with an error. After each error within
%                           the try-catch block, the function
%                           cleanupMatlabPrefs is called. This is thought
%                           to correct a stochastic error that can occur in
%                           parpool jobs and results in the corruption of
%                           the matlab preference file, which is then
%                           deleted.
%   'maxAttempts'         - The number of times that a given stage will be 
%                           re-tried in the event of an error.
%   'makeFitVideoByName'  - A cell array of stages for which a fit
%                           video will be produced following completion of
%                           the stage. The video file name is the run name,
%                           followed by "_fitStageX.avi" where X is the idx
%                           of the function call list.
%   'makeFitVideoByNumber' - A cell array of stages for which a fit
%                           video will be produced following completion of
%                           the stage. The video file name is the run name,
%                           followed by "_fitStageX.avi" where X is the idx
%                           of the function call list.
%
% Output:
%   None. The routine saves files but does not return variables.
%


%% Parse input and define variables
p = inputParser; p.KeepUnmatched = true;

% required input
p.addRequired('pathParams',@isstruct);

% optional input
p.addParameter('lastStageByName', {}, @iscell);
p.addParameter('lastStageByNumber', [], @isnumeric);
p.addParameter('skipStageByName', {}, @iscell);
p.addParameter('skipStageByNumber',[], @isnumeric);
p.addParameter('displayAvailableStages', false, @islogical)
p.addParameter('rawVideoSuffix', {'_raw.mov' '.mov'}, @iscell);
p.addParameter('videoTypeChoice', 'LiveTrackWithVTOP_eye', @ischar);
p.addParameter('customFunCalls', {}, @iscell);
p.addParameter('customSceneGeometryFile', [], @(x)(isempty(x) | iscell(x) | ischar(x)));
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

% Test if there are instructions to skip the deinterlaceVideo stage. If
% not, try to identify the raw video, which may be in one of several video
% formats or with different suffix forms.
if ~any(strcmp(p.Results.skipStageByName,'deinterlaceVideo')) && ~any(p.Results.skipStageByNumber == 1) 
    % Create a cell array of candidate raw video names with the runName and
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
if isfield(pathParams,'grayVideoName')
    grayVideoName = pathParams.grayVideoName;
else
    grayVideoName = fullfile(pathParams.dataOutputDirFull, [pathParams.runName '_gray.avi']);
end
timebaseFileName = fullfile(pathParams.dataOutputDirFull, [pathParams.runName '_timebase.mat']);
glintFileName = fullfile(pathParams.dataOutputDirFull, [pathParams.runName '_glint.mat']);
perimeterFileName = fullfile(pathParams.dataOutputDirFull, [pathParams.runName '_perimeter.mat']);
controlFileName = fullfile(pathParams.dataOutputDirFull, [pathParams.runName '_controlFile.csv']);
correctedPerimeterFileName = fullfile(pathParams.dataOutputDirFull, [pathParams.runName '_correctedPerimeter.mat']);
pupilFileName = fullfile(pathParams.dataOutputDirFull, [pathParams.runName '_pupil.mat']);

% The sceneGeometryFileNameOutput is the name of the file created by the
% estimateSceneParams routine.
sceneGeometryFileNameOutput = fullfile(pathParams.dataOutputDirFull, [pathParams.runName '_sceneGeometry.mat']);

% The sceneGeometryFileNameInput is the name of the file passed to the
% fitPupilPerimeter, smoothPupilRadius, and makeFitVideo routines. A custom
% value can be passed. The custom value will usually identify a different
% acqusition during which the head was in the same positiion, and for which
% sceneGeometry was well estimated.
if ~isempty(p.Results.customSceneGeometryFile)
    if iscell(p.Results.customSceneGeometryFile)
        sceneGeometryFileNameInput = p.Results.customSceneGeometryFile{1};
    else
        sceneGeometryFileNameInput = p.Results.customSceneGeometryFile;
    end
else
    sceneGeometryFileNameInput = fullfile(pathParams.dataOutputDirFull, [pathParams.runName '_sceneGeometry.mat']);
end
sceneDiagnosticPlotFileName = fullfile(pathParams.dataOutputDirFull, [pathParams.runName '_sceneDiagnosticPlot.pdf']);
finalFitVideoName = fullfile(pathParams.dataOutputDirFull, [pathParams.runName '_finalFit.avi']);
eyeModelVideoName = fullfile(pathParams.dataOutputDirFull, [pathParams.runName '_eyeModel.avi']);


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
            ' ''pupilFileName'', pupilFileName, ''fitLabel'', ''initial'',' ...
            ' varargin{:});']...
            };
    case 'LiveTrackWithVTOP_eye'
        funCalls = {...
            'deinterlaceVideo(rawVideoName, grayVideoName, ''timebaseFileName'',timebaseFileName,varargin{:});'...
            'findGlint(grayVideoName, glintFileName, varargin{:});'...
            'findPupilPerimeter(grayVideoName, perimeterFileName, varargin{:});'...
            'makeControlFile(controlFileName, perimeterFileName, glintFileName, varargin{:});' ...
            'applyControlFile(perimeterFileName,controlFileName,correctedPerimeterFileName, varargin{:});' ...
            ['fitPupilPerimeter(correctedPerimeterFileName, pupilFileName,' ...
            '''nSplits'', 0, varargin{:});']...
            ['estimateSceneParams(pupilFileName, sceneGeometryFileNameOutput,' ...
            '''sceneDiagnosticPlotFileName'', sceneDiagnosticPlotFileName, varargin{:});']...
            ['fitPupilPerimeter(correctedPerimeterFileName, pupilFileName,' ...
            '''sceneGeometryFileName'', sceneGeometryFileNameInput, varargin{:});']...
            'smoothPupilRadius(correctedPerimeterFileName, pupilFileName, sceneGeometryFileNameInput, varargin{:});'...
            ['makeFitVideo(grayVideoName, finalFitVideoName,' ...
            '''glintFileName'', glintFileName, ''perimeterFileName'', correctedPerimeterFileName,'...
            '''pupilFileName'', pupilFileName, ''sceneGeometryFileName'', sceneGeometryFileNameInput,' ...
            '''modelEyeAlpha'', 0.25,' ...
            '''controlFileName'',controlFileName,varargin{:});']...
            'makeEyeModelVideo(eyeModelVideoName, pupilFileName, sceneGeometryFileNameInput, varargin{:});'...
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
        
    end % if we aren't skipping this stage by name or number
    
    % Check if we should make a fit video for this stage. Note that it is
    % possible that we would make a fit video even though we skipped the
    % stage
    if any(strcmp(p.Results.makeFitVideoByName,funNames{ff})) || any(p.Results.makeFitVideoByNumber == ff)
        makeFitVideoForThisStage(pathParams, sceneGeometryFileNameInput, funNames, ff, varargin{:});
    end
    
    % clear all files (hopefully prevents 'too many files open' error)
    fclose all;
    
    % Check to see if we have just completed the last stage, in which case
    % we should break out of the loop of function calls
    if strcmp(p.Results.lastStageByName,funNames{ff})
        break
    end
    if p.Results.lastStageByNumber==ff
        break
    end
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


function makeFitVideoForThisStage(pathParams, sceneGeometryFileNameInput, funNames, ff, varargin)

% Define the fitVideo output name
fitVideoFileName = fullfile(pathParams.dataOutputDirFull, [pathParams.runName '_fitStage' num2str(ff) '.avi']);

% Assemble the entire list of potential files to include in the video
if isfield(pathParams,'grayVideoName')
    grayVideoName = pathParams.grayVideoName;
else
    grayVideoName = fullfile(pathParams.dataOutputDirFull, [pathParams.runName '_gray.avi']);
end
glintFileName = fullfile(pathParams.dataOutputDirFull, [pathParams.runName '_glint.mat']);
initialPerimeterFileName = fullfile(pathParams.dataOutputDirFull, [pathParams.runName '_perimeter.mat']);
controlFileName = fullfile(pathParams.dataOutputDirFull, [pathParams.runName '_controlFile.csv']);
correctedPerimeterFileName = fullfile(pathParams.dataOutputDirFull, [pathParams.runName '_correctedPerimeter.mat']);
pupilFileName = fullfile(pathParams.dataOutputDirFull, [pathParams.runName '_pupil.mat']);

% Depending upon which stage just completed, set to empty any files
% that do not yet exist, and select the appropriate perimeter file and
% pupil ellipse fit to display.
switch funNames{ff}
    case 'deinterlaceVideo'
        glintFileName = []; perimeterFileName=[]; controlFileName=[]; pupilFileName=[]; sceneGeometryFileNameInput=[];
    case 'resizeAndCropVideo'
        glintFileName = []; perimeterFileName=[]; controlFileName=[]; pupilFileName=[]; sceneGeometryFileNameInput=[];
    case 'findGlint'
        perimeterFileName=[]; controlFileName=[]; pupilFileName=[]; sceneGeometryFileNameInput=[];
    case 'findPupilPerimeter'
        if ~exist(glintFileName,'file')
            glintFileName = [];
        end
        perimeterFileName=initialPerimeterFileName;
        controlFileName=[]; pupilFileName=[]; sceneGeometryFileNameInput=[];
    case 'makeControlFile'
        if ~exist(glintFileName,'file')
            glintFileName = [];
        end
        perimeterFileName=initialPerimeterFileName;
        controlFileName=[]; pupilFileName=[]; sceneGeometryFileNameInput=[];
    case 'applyControlFile'
        if ~exist(glintFileName,'file')
            glintFileName = [];
        end
        perimeterFileName=correctedPerimeterFileName;
        pupilFileName=[]; sceneGeometryFileNameInput=[];
    case 'fitPupilPerimeter'
        if ~exist(glintFileName,'file')
            glintFileName = [];
        end
        % If a control file has not been created, set this to empty
        controlFunCallIdx=find(strcmp(funNames,'makeControlFile'));
        if isempty(controlFunCallIdx)
            controlFileName = [];
        else
            if controlFunCallIdx > ff
                controlFileName = [];
            end
        end
        % Show the initial perimeter, unless it has been corrected
        perimeterFileName=initialPerimeterFileName;
        correctFunCallIdx=find(strcmp(funNames,'applyControlFile'));
        if ~isempty(correctFunCallIdx)
            if correctFunCallIdx < ff
                perimeterFileName=correctedPerimeterFileName;
            end
        end
        % If the sceneGeometry has been determined by this point, we can
        % plot the sceneConstrained ellipse fit, otherwise plot the
        % unconstrained
        varargin={varargin{:}, 'fitLabel', 'initial'};
        sceneFunCallIdx=find(strcmp(funNames,'estimateSceneParams'));
        if ~isempty(sceneFunCallIdx)
            if sceneFunCallIdx < ff
                varargin={varargin{:}, 'fitLabel', 'sceneConstrained'};
            else
                sceneGeometryFileNameInput=[];
            end
        else
            sceneGeometryFileNameInput=[];
        end
    case 'estimateSceneParams'
        if ~exist(glintFileName,'file')
            glintFileName = [];
        end
        perimeterFileName=correctedPerimeterFileName;
        varargin={varargin{:}, 'fitLabel', 'initial'};
    case 'smoothPupilRadius'
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
    'sceneGeometryFileName', sceneGeometryFileNameInput, varargin{:});

end % makeFitVideoForThisStage




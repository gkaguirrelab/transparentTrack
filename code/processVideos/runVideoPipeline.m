function runVideoPipeline( pathParams, varargin )
% runVideoPipeline( pathParams, varargin )
%
% this is the standard processing pipeline for eye tracking videos.
%
% The pipeline consists of the following stages:
%   deinterlaceVideo
%   resizeAndCropVideo
%   findGlint
%   findPupilPerimeter
%   makeControlFile
%   applyControlFile
%   fitPupilPerimeter
%   smoothPupilParameters
%   fitIrisPerimeter
%   makeFitVideo
%
% The user can stop the execution after any of the stages with the optional
% param 'lastStage', or skip any number of stages listing them in a cell
% under the optional param 'skipStage'. Every stage, however, requires the
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
%   skipStage - a cell array of function calls to be skipped during
%      execution of the pipeline.
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
%


%% Parse input and define variables
p = inputParser; p.KeepUnmatched = true;

% required input
p.addRequired('pathParams',@isstruct);

% optional input
p.addParameter('lastStage', 'makePupilFitVideo', @ischar);
p.addParameter('skipStage', {}, @iscell);
p.addParameter('rawVideoSuffix', {'_raw.mov' '.mov'}, @iscell);
p.addParameter('videoTypeChoice', 'LiveTrackWithVTOP_eyeNoIris', @ischar);
p.addParameter('customFunCalls', {}, @iscell);
p.addParameter('catchErrors', true, @islogical);

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
if ~exist(pathParams.controlFileDirFull,'dir')
    mkdir(pathParams.controlFileDirFull)
end


%% Define input and output filenames

if ~any(strcmp(p.Results.skipStage,'convertRawToGray'))
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
grayVideoName = fullfile(pathParams.dataOutputDirFull, [pathParams.runName '_gray.avi']);
glintFileName = fullfile(pathParams.dataOutputDirFull, [pathParams.runName '_glint.mat']);
perimeterFileName = fullfile(pathParams.dataOutputDirFull, [pathParams.runName '_perimeter.mat']);
controlFileName = fullfile(pathParams.controlFileDirFull, [pathParams.runName '_controlFile.csv']);
correctedPerimeterFileName = fullfile(pathParams.dataOutputDirFull, [pathParams.runName '_correctedPerimeter.mat']);
pupilFileName = fullfile(pathParams.dataOutputDirFull, [pathParams.runName '_pupil.mat']);
sceneGeometryFileName = fullfile(pathParams.dataOutputDirFull, [pathParams.runName '_sceneGeometry.mat']);
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
            ' ''nsplits'', 1, varargin{:});']...
            ['smoothPupilParameters(perimeterFileName, pupilFileName,' ...
            ' ''exponentialTauParams'', [100, 100, 100, 100, 100], varargin{:});']...
            ['makeFitVideo(grayVideoName, finalFitVideoName,' ...
            ' ''perimeterFileName'', perimeterFileName,''perimeterColor'',''r'','...
            ' ''pupilFileName'', pupilFileName, ''whichFieldToPlot'', ''pPosteriorMeanTransparent'',' ...
            ' varargin{:});']...
            };
    case 'LiveTrackWithVTOP_eye'
        funCalls = {...
            'deinterlaceVideo(rawVideoName,grayVideoName, varargin{:});'...
            'findGlint(grayVideoName, glintFileName, varargin{:});'...
            'findPupilPerimeter(grayVideoName, perimeterFileName, varargin{:});'...
            'makeControlFile(controlFileName, perimeterFileName, glintFileName, varargin{:});' ...
            'applyControlFile(perimeterFileName,controlFileName,correctedPerimeterFileName, varargin{:});' ...
            'fitPupilPerimeter(correctedPerimeterFileName, pupilFileName, varargin{:});'...
            'estimateSceneGeometry(pupilFileName, sceneGeometryFileName, varargin{:});'...
            ['makeControlFile(controlFileName, perimeterFileName, glintFileName,' ...
            '''sceneGeometryFileName'', sceneGeometryFileName, ''overwriteControlFile'', true, varargin{:});']...
            ['fitPupilPerimeter(correctedPerimeterFileName, pupilFileName,' ...
            '''sceneGeometryFileName'', sceneGeometryFileName, varargin{:});']...
            'smoothPupilParameters(correctedPerimeterFileName, pupilFileName, varargin{:});'...
            'fitIrisPerimeter(grayVideoName, perimeterFileName, pupilFileName, irisFileName, varargin{:});' ...
            ['makeFitVideo(grayVideoName, finalFitVideoName,' ...
            '''glintFileName'', glintFileName, ''perimeterFileName'', correctedPerimeterFileName,'...
            '''pupilFileName'', pupilFileName, ''whichFieldToPlot'', ''pPosteriorMeanTransparent'',' ...
            '''irisFileName'', irisFileName,' ...
            '''controlFileName'',controlFileName,varargin{:});']...
            };
        case 'LiveTrackWithVTOP_eyeNoIris'
        funCalls = {...
            'deinterlaceVideo(rawVideoName,grayVideoName, varargin{:});'...
            'findGlint(grayVideoName, glintFileName, varargin{:});'...
            'findPupilPerimeter(grayVideoName, perimeterFileName, varargin{:});'...
            'makeControlFile(controlFileName, perimeterFileName, glintFileName, varargin{:});' ...
            'applyControlFile(perimeterFileName,controlFileName,correctedPerimeterFileName, varargin{:});' ...
            'fitPupilPerimeter(correctedPerimeterFileName, pupilFileName, varargin{:});'...
            'estimateSceneGeometry(pupilFileName, sceneGeometryFileName, varargin{:});'...
            ['makeControlFile(controlFileName, perimeterFileName, glintFileName,' ...
            '''sceneGeometryFileName'', sceneGeometryFileName, ''overwriteControlFile'', true, varargin{:});']...
            ['fitPupilPerimeter(correctedPerimeterFileName, pupilFileName,' ...
            '''sceneGeometryFileName'', sceneGeometryFileName, varargin{:});']...
            'smoothPupilParameters(correctedPerimeterFileName, pupilFileName, varargin{:});'...
            ['makeFitVideo(grayVideoName, finalFitVideoName,' ...
            '''glintFileName'', glintFileName, ''perimeterFileName'', correctedPerimeterFileName,'...
            '''pupilFileName'', pupilFileName, ''whichFieldToPlot'', ''pPosteriorMeanTransparent'',' ...
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

% Loop through the function calls
for ff = 1:length(funCalls)
    if ~any(strcmp(p.Results.skipStage,funNames{ff}))

        % check if we are going to catch or throw errors
        if p.Results.catchErrors
            % intialize while control
            success = 0;
            attempts = 1; % we will attempt to execute the instruction 3 times. If it does not work, the code will eventually break.
            while ~success
                try
                    eval(funCalls{ff});
                    success = 1;
                catch ME
                    warning ('There has been an error during execution. Cleaning matlabprefs.mat and trying again - attempt %d.',attempts)
                    cleanupMatlabPrefs;
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
        else % we will throw errors
            eval(funCalls{ff});
        end % if catchErrors

        % clear all files (hopefully prevents 'too many files open' error)
        fclose all ;
        if strcmp(p.Results.lastStage,funNames{ff})
            return
        end
    end % if we aren't skipping this stage
end % loop over function calls

end % main function



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

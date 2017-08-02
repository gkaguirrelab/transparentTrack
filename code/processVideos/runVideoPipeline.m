function runVideoPipeline( pathParams, varargin )
% runVideoPipeline( pathParams, varargin )
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
%   fitIrisCircleAndMask
%   makePupilFitVideo
%
% The user can stop the execution after any of the stages with the optional
% param 'lastStage', or skip any amount of stages listing them in a cell
% under the optional param 'skipStage'. Every stage however requires the
% existence of the output from the preceeding ones to be correctly
% executed.
%
% INPUT
%   pathParams - DESCRIBE THIS STRUCTURE AND ITS FIELDS HERE
%
% OUTPUT
%   None. The routine saves files but does not return variables.
%
% OPTIONS
%   rawVideoSuffix - cell array of strings that contain possible suffixes
%       of raw video files to be processed. The routine will search for
%       raw videos sequentially in the cell array until it finds a match
%

%% Parse input and define variables
p = inputParser; p.KeepUnmatched = true;

% required input
p.addRequired('pathParams',@isstruct);

% optional input
p.addParameter('lastStage', 'makePupilFitVideo', @ischar);
p.addParameter('rawVideoSuffix', {'_raw.mov' '.mov'}, @iscell);
p.addParameter('videoTypeChoice', 'LiveTrackWithVTOP_eye', @ischar);
p.addParameter('customFunCalls', {}, @iscell);
p.addParameter('skipStage', {}, @iscell);

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

grayVideoName = fullfile(pathParams.dataOutputDirFull, [pathParams.runName '_gray.avi']);
glintFileName = fullfile(pathParams.dataOutputDirFull, [pathParams.runName '_glint.mat']);
perimeterFileName = fullfile(pathParams.dataOutputDirFull, [pathParams.runName '_perimeter.mat']);
controlFileName = fullfile(pathParams.controlFileDirFull, [pathParams.runName '_controlFile.csv']);
correctedPerimeterFileName = fullfile(pathParams.dataOutputDirFull, [pathParams.runName '_correctedPerimeter.mat']);
pupilFileName = fullfile(pathParams.dataOutputDirFull, [pathParams.runName '_pupil.mat']);
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
            'convertRawToGray(rawVideoName,grayVideoName, varargin{:});'...
            'findPupilPerimeter(grayVideoName, perimeterFileName, ''pupilOnly'', true, varargin{:});'...
            'fitPupilPerimeter(perimeterFileName, pupilFileName, varargin{:});'...
            ['makeFitVideo(grayVideoName, finalFitVideoName,' ...
            ' ''perimeterFileName'', perimeterFileName,''perimeterColor'',''r'','...
            ' ''pupilFileName'', pupilFileName, ''whichFieldToPlot'', ''pPosteriorMeanTransparent'',' ...
            ' varargin{:});']...
            };
    case 'LiveTrackWithVTOP_eye'
        funCalls = {...
            'convertRawToGray(rawVideoName,grayVideoName, varargin{:});'...
            'trackGlint(grayVideoName, glintFileName, varargin{:});'...
            'findPupilPerimeter(grayVideoName, perimeterFileName, varargin{:});'...
            'makeControlFile(controlFileName, perimeterFileName, glintFileName, varargin{:});' ...
            'applyControlFile(perimeterFileName,controlFileName,correctedPerimeterFileName, varargin{:});' ...
            'fitPupilPerimeter(correctedPerimeterFileName, pupilFileName, varargin{:});'...
            'fitIrisPerimeter(grayVideoName, perimeterFileName, pupilFileName, irisFileName, varargin{:});' ...
            ['makeFitVideo(grayVideoName, finalFitVideoName,' ...
            '''glintFileName'', glintFileName, ''perimeterFileName'', correctedPerimeterFileName,'...
            '''pupilFileName'', pupilFileName, ''whichFieldToPlot'', ''pPosteriorMeanTransparent'',' ...
            '''irisFileName'', irisFileName,' ...
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
        if strcmp(p.Results.lastStage,funNames{ff})
            return
        end
    end
end

end % main function



function cleanupMatlabPrefs    
% this small function closes any open parallel pool and removes the
% matlabprefs.mat file in case it is corrupt.
% This file can get corrupt while using parallel pools because all of the
% pool worker try to read/write it at the same time. MathWorks has been
% contacted about this issue.

% close parallel pool if any is open
poolobj = gcp('nocreate');
if ~isempty (poolobj)
    delete(poolobj);
end

% clean up the corrupt matlabprefs file
prefFile = fullfile(prefdir,'matlabprefs.mat');
system(['rm -rf "' prefFile '"'])
end % cleanupMatlabPrefs

function makePreliminaryControlFile(controlFileName, perimeterFileName, glintFileName, varargin)
% function makePreliminaryControlFile(controlFileName, perimeterFileName, glintFileName, varargin)
%
% The routine creates and saves a "control file", which is a text (csv)
% file that instruct subsequent routines as to how a perimeterFile may be
% cleaned up to remove blinks, eyelid intrusions, or force a particular
% set of ellipse parameters.
%
% The control file created here is an initial guess at these instructions.
% It includes two primary routnes:
%
% 1) Blink detection. Frames that lack a glint entirely, or in which the
% glint is at a non-plausible location when compared to the average glint
% location, are marked as blinks. Frames that are adjacent to blinks may
% also be excluded.
% 2) A more time-consuming step examines if removing a portion of the
% perimeter of the pupil boundary produces an improvement in an initial
% ellipse fit. This routine operates in a parfor loop.
%
%
% FORMAT OF CONTROL FILE:
%	Each line of the control file contains "instruction" of the form:
%       FRAME NUMBER, INSTRUCTION TYPE, INSTRUCTION PARAMS (variable #)
%
%	where:
%       FRAME NUMBER: frame on which to apply the instruction.
%       INSTRUCTION TYPE: what to do on the frame.
%       INSTRUCTION PARAMS: variable number of params necessary to execute
%           the instruction.
%
%   This routine generates instructions of the type "blink" and "cut".
%   Other available instructions are "ellipse", with the parameters of a
%   forced ellipse, and "%" which marks a non-executed comment
%
% Input (required)
%	controlFileName - full path to the control file (including csv extension)
%   perimeterFileName -
%   glintFileName - 
%
% Options (analysis)
%   glintDisplaceSTD - Glint displacemenet more than this number of
%       standard deviations away from the mean glint position will be
%       judged to be a reflection off the eyelid of the subject and thus
%       part of a blink.
%	extendBlinkWindow - a two element vector that defines the number of
%       additional frames flagged as a blink before and after a continuous
%       block blinks.
%   cutErrorThreshold - the distance error tolerated before attempting to
%       cut
%   ellipseTransparentLB/UB - the lower and upper bounds of the constrained
%      ellipse fit that is used to judge the quality of different cuts.
%
% Optional key/value pairs (verbosity and I/O)
%  'verbosity' - level of verbosity. [none, full]
%  'ellipseFitDataFileName': full path to the .mat file in which to save
%     pupil tracking information.
%
% Optional key/value pairs (flow control)
%  'useParallel' - If set to true, use the Matlab parallel pool for the
%    initial ellipse fitting.
%  'nWorkers' - Specify the number of workers in the parallel pool. If
%    undefined the default number will be used.
%  'tbtbProjectName' - The workers in the parallel pool are configured by
%    issuing a tbUseProject command for the project specified here.
%
% Optional key/value pairs (Environment parameters)
%  'tbSnapshot' - This should contain the output of the tbDeploymentSnapshot
%    performed upon the result of the tbUse command. This documents the
%    state of the system at the time of analysis.
%  'timestamp' - AUTOMATIC - The current time and date
%  'username' - AUTOMATIC - The user
%  'hostname' - AUTOMATIC - The host


%% Parse input and define variables
p = inputParser;

% required input
p.addRequired('controlFileName',@isstr);
p.addRequired('perimeterFileName',@isstr);
p.addRequired('glintFileName',@isstr);

% Optional analysis params
p.addParameter('glintDisplaceSTD', 2, @isnumeric);
p.addParameter('extendBlinkWindow', [2,2], @isnumeric);
p.addParameter('cutErrorThreshold', 10, @isnumeric);
p.addParameter('ellipseTransparentLB',[0, 0, 1000, 0, -0.5*pi],@isnumeric);
p.addParameter('ellipseTransparentUB',[240,320,10000,0.417, 0.5*pi],@isnumeric);

% Optional display params
p.addParameter('verbosity','none',@ischar);

% Optional flow control params
p.addParameter('useParallel',true,@islogical);
p.addParameter('nWorkers',[],@(x)(isempty(x) | isnumeric(x)));
p.addParameter('tbtbRepoName','LiveTrackAnalysisToolbox',@ischar);

% Environment parameters
p.addParameter('tbSnapshot',[],@(x)(isempty(x) | isstruct(x)));
p.addParameter('timestamp',char(datetime('now')),@ischar);
p.addParameter('username',char(java.lang.System.getProperty('user.name')),@ischar);
p.addParameter('hostname',char(java.net.InetAddress.getLocalHost.getHostName),@ischar);

% parse
p.parse(controlFileName, perimeterFileName, glintFileName, varargin{:})


%% Check that this control file does not exist
% We decline to over-write an existing control file, as it may contain
% instructions lovingly crafted by a human
if exist(controlFileName, 'file') == 2
    error(['The control file ' controlFileName ' exists already and cannot be over-written']);
end


%% Set up the parallel pool
if p.Results.useParallel
    if strcmp(p.Results.verbosity,'full')
        tic
        fprintf(['Opening parallel pool. Started ' char(datetime('now')) '\n']);
    end
    if isempty(p.Results.nWorkers)
        parpool;
    else
        parpool(p.Results.nWorkers);
    end
    poolObj = gcp;
    if isempty(poolObj)
        nWorkers=0;
    else
        nWorkers = poolObj.NumWorkers;
        % Use TbTb to configure the workers.
        if ~isempty(p.Results.tbtbRepoName)
            spmd
                tbUse(p.Results.tbtbRepoName,'reset','full','verbose',false,'online',false);
            end
            if strcmp(p.Results.verbosity,'full')
                fprintf('CAUTION: Any TbTb messages from the workers will not be shown.\n');
            end
        end
    end
    if strcmp(p.Results.verbosity,'full')
        toc
        fprintf('\n');
    end
else
    nWorkers=0;
end


%% Perform blink detection
dataLoad = load(glintFileName);
glintData = dataLoad.glintData;
clear dataLoad

% locate all nans
blinks = find(isnan(glintData.X));

% locate all points further away than glintDisplaceSTD from the glint mean
farX = find(glintData.X > nanmean(glintData.X) + (p.Results.glintDisplaceSTD * nanstd(glintData.X)) | ...
    glintData.X < nanmean(glintData.X) - (p.Results.glintDisplaceSTD * nanstd(glintData.X)));
farY = find(glintData.Y > nanmean(glintData.Y) + (p.Results.glintDisplaceSTD * nanstd(glintData.Y)) | ...
    glintData.Y < nanmean(glintData.Y) - (p.Results.glintDisplaceSTD * nanstd(glintData.Y)));
farGlints = union(farX, farY);

% combine to obtain all blinks
blinkFrames = union(blinks, farGlints);

% ADD CODE HERE TO EXTEND NANS ACCORDING TO THE extendBlinkWindow param


%% Guess pupil cuts

% Load the pupil perimeter data. It will be a structure variable
% "perimeter", with the fields .data and .meta
dataLoad=load(perimeterFileName);
perimeter=dataLoad.perimeter;
clear dataLoad
nFrames=size(perimeter.data,3);

% Intialize some variables
frameRadii=nan(nFrames,1);
frameThetas=nan(nFrames,1);

% alert the user
if strcmp(p.Results.verbosity,'full')
    tic
    fprintf(['Guessing pupil cuts. Started ' char(datetime('now')) '\n']);
    fprintf('| 0                      50                   100%% |\n');
    fprintf('.\n');
end

% Loop through the video frames
parfor (ii = 1:nFrames, nWorkers)
    
    % Update progress
    if strcmp(p.Results.verbosity,'full') && mod(ii,round(nFrames/50))==0
        fprintf('\b.\n');
    end
    
    % get the data frame
    binP = squeeze(perimeter.data(:,:,ii));
    [Yp, Xp] = ind2sub(size(binP),find(binP));
    
    % proceed if the frame is not empty and has not been tagged as a blink
    if ~ismember(ii,blinkFrames) && ~isempty(Xp)
        
        % fit an ellipse to the full perimeter using the constrainedEllipseFit
        [~, ~, originalFittingError] = constrainedEllipseFit(Xp, Yp, ...
            p.Results.ellipseTransparentLB, ...
            p.Results.ellipseTransparentUB, []);
        
        % if the fitting error is above the threshold, search over cuts
        if originalFittingError > p.Results.cutErrorThreshold
            candidateRadii=0;
            candidateThetas=pi/2:pi/16:pi;
            [gridSearchRadii,gridSearchThetas] = ndgrid(candidateRadii,candidateThetas);
            myCutOptim = @(params) calcErrorForACut(binP, params(1), params(2), p.Results.ellipseTransparentLB, p.Results.ellipseTransparentUB);
            gridSearchResults=arrayfun(@(k1,k2) myCutOptim([k1,k2]),gridSearchRadii,gridSearchThetas);
            
            [row,col] = find(gridSearchResults==min(min(gridSearchResults)));
            
            if ~isempty(row)
                frameRadii(ii)=candidateRadii(row(1));
                frameThetas(ii)=candidateThetas(col(1));
            end
            
        end % search over cuts
        
    end % not an empty frame
    
end % parloop over frames

% report completion of preliminary control file generation
if strcmp(p.Results.verbosity,'full')
    toc
    fprintf('\n');
end

%% Write out control file instructions
fid = fopen(controlFileName,'a');
% write out blinks
if ~isempty(blinkFrames)
for bb = 1 : length(blinkFrames)
    instruction = [num2str(blinkFrames(bb)) ',' 'blink'];
    fprintf(fid,'%s\n',instruction);
    clear instruction
end
end
% write out cuts
cutFrames=find(~isnan(frameThetas));
if ~isempty(cutFrames)
    for cc = 1 : length(cutFrames)
        frameIdx=cutFrames(cc);
        instruction = [num2str(frameIdx) ',' 'cut' ',' num2str(frameRadii(frameIdx)) ',' num2str(frameThetas(frameIdx))];
        fprintf(fid,'%s\n',instruction);
        clear instruction
    end
end
% finish and close the file
instruction = ['%' ',' '%' ',' 'end of automatic instructions'];
fprintf(fid,'%s\n',instruction);
fclose(fid);

%% Delete the parallel pool
if strcmp(p.Results.verbosity,'full')
    tic
    fprintf(['Closing parallel pool. Started ' char(datetime('now')) '\n']);
end
if p.Results.useParallel
    poolObj = gcp;
    if ~isempty(poolObj)
        delete(poolObj);
    end
end
if strcmp(p.Results.verbosity,'full')
    toc
    fprintf('\n');
end


end % function

function distanceError = calcErrorForACut(theFrame, radiusThresh, theta, lb, ub)
[binPcut] = cutPupil (theFrame, radiusThresh, theta);
[Yp, Xp] = ind2sub(size(binPcut),find(binPcut));
[~, ~, distanceError] = constrainedEllipseFit(Xp, Yp, lb, ub, []);
end

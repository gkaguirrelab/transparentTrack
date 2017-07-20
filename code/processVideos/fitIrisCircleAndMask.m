function fitIrisCircleAndMask(grayVideoName, perimeterFileName, pupilFitFileName, irisFitFileName, varargin)
% function fitIrisAndPalpebralFissure(grayVideoName, perimeterFileName, pupilFitFileName, irisFitFileName, palpebralFissureFileName, varargin)
%
% This function fits a circle to the outer border of the iris, and creates
% a palpebral fissure mask for each frame. The core of the routine makes
% use of the IrisSeg toolbox (https://github.com/cdac-cvml/IrisSeg).
%
% Output
%   irisFitData - structure with a 'data' field containing the parameters
%       of ellipse fits to the border of the iris, and a meta field with
%       analysis and environment params.
%
% Input (required)
%	grayVideoName - full path to  the gray video to track
%	ellipseFitDataFileName - full path to the .mat file that contains fits
%       to the border of the pupil
%
% Options (analysis)
% 	gammaCorrection - gamma correction to be applied to the video frames
%       (default 1, typical range [0.5 1.8])
%
% Options (verbosity and display)
%   verbosity - controls console status updates
%
% Optional key/value pairs (flow control)
%
%  'nFrames' - analyze fewer than the total number of frames.
%  'useParallel' - If set to true, use the Matlab parallel pool for the
%    initial ellipse fitting.
%  'nWorkers' - Specify the number of workers in the parallel pool. If
%    undefined the default number will be used.
%  'tbtbProjectName' - The workers in the parallel pool are configured by
%    issuing a tbUseProject command for the project specified here.
%  'developmentMode' - If set to true, the routine attempts to load a
%    pre-existing set of initial ellipse measures (and SDs upon those
%    params), rather than re-computing these. This allows more rapid
%    exploration of parameter settigns that guide the Bayesian smoothing.
%
% Optional key/value pairs (Environment parameters)
%  'tbSnapshot' - This should contain the output of the tbDeploymentSnapshot
%    performed upon the result of the tbUse command. This documents the
%    state of the system at the time of analysis.
%  'timestamp' - AUTOMATIC - The current time and date
%  'username' - AUTOMATIC - The user
%  'hostname' - AUTOMATIC - The host
%
%


%% Parse input and define variables
p = inputParser; p.KeepUnmatched = true;

% required input
p.addRequired('grayVideoName',@isstr);
p.addRequired('perimeterFileName',@isstr);
p.addRequired('pupilFitFileName',@isstr);
p.addRequired('irisFitFileName',@isstr);

% Optional display params
p.addParameter('verbosity','none',@ischar);

% Optional analysis params
p.addParameter('gammaCorrection', 1, @isnumeric);

% Optional flow control params
p.addParameter('nFrames',Inf,@isnumeric);
p.addParameter('useParallel',false,@islogical);
p.addParameter('nWorkers',[],@(x)(isempty(x) | isnumeric(x)));
p.addParameter('tbtbRepoName','LiveTrackAnalysisToolbox',@ischar);

% Environment parameters
p.addParameter('tbSnapshot',[],@(x)(isempty(x) | isstruct(x)));
p.addParameter('timestamp',char(datetime('now')),@ischar);
p.addParameter('username',char(java.lang.System.getProperty('user.name')),@ischar);
p.addParameter('hostname',char(java.net.InetAddress.getLocalHost.getHostName),@ischar);

% parse
p.parse(grayVideoName, perimeterFileName, pupilFitFileName, irisFitFileName, varargin{:})


%% Read files into memory

% load gray video
videoInObj = VideoReader(grayVideoName);
% get number of frames
if p.Results.nFrames == Inf
    nFrames = floor(videoInObj.Duration*videoInObj.FrameRate);
else
    nFrames = p.Results.nFrames;
end
% get video dimensions
videoSizeX = videoInObj.Width;
videoSizeY = videoInObj.Height;
% initialize variable to hold the perimeter data
grayVideo = zeros(videoSizeY,videoSizeX,nFrames,'uint8');
% read the video into memory, adjusting gamma and local contrast
for ii = 1:min([floor(videoInObj.Duration*videoInObj.FrameRate) p.Results.nFrames])
    thisFrame = readFrame(videoInObj);
    thisFrame = imadjust(thisFrame,[],[],p.Results.gammaCorrection);
    grayVideo(:,:,ii) = rgb2gray (thisFrame);
end
% close the video object
clear videoInObj

% Read in the ellipseFitData
dataLoad = load(pupilFitFileName);
pupilFitData = dataLoad.ellipseFitData;
clear dataLoad
pupilFitParams = pupilFitData.pPosteriorMeanTransparent;

% Load the pupil perimeter data. It will be a structure variable
% "perimeter", with the fields .data and .meta
dataLoad=load(perimeterFileName);
pupilPerimeter=dataLoad.perimeter;
clear dataLoad


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


%% Fit iris perimeter

% Detect display mode

% alert the user
if strcmp(p.Results.verbosity,'full')
    tic
    fprintf(['Initial iris width detection. Started ' char(datetime('now')) '\n']);
    fprintf('| 0                      50                   100%% |\n');
    fprintf('.\n');
end

% initialize variables to hold the results
irisInitialWidth = nan(nFrames,1);

% Initial loop through gray frames to obtain an estimate of iris width
parfor (ii = 1:nFrames, nWorkers)
    
    % increment the progress bar
    if strcmp(p.Results.verbosity,'full') && mod(ii,round(nFrames/50))==0
        fprintf('\b.\n');
    end
    
    % diagnostic try/catch
    try
        
        % check if there is a defined pupil fit. If so, proceed
        if ~isnan(pupilFitParams(ii,1))
            
            % get the frame
            thisFrame = squeeze(grayVideo(:,:,ii));
            
            % Calculate the pupil height from the corrected pupil perimeter
            pupil_height = min(find(~max(squeeze(pupilPerimeter.data(:,:,ii))')==0));
            
            % To maintain transparency with regard to the IrisSeg toolbox, we
            % transfer values to variable names used by those routines
            scale= 1;
            pCentreX=pupilFitParams(ii,1);
            pCentreY=pupilFitParams(ii,2);
            pRadius=sqrt(pupilFitParams(ii,3)/pi);
            N = videoSizeY;
            M = videoSizeX;
            
            %% Code pulled from irisseg_main.m, part of the IrisSeg toolbox
            % https://github.com/cdac-cvml/IrisSeg
            
            AngRadius1 = pRadius * 8;
            if ( (AngRadius1 > pCentreX) || (AngRadius1 > pCentreY) || (AngRadius1 > (M-pCentreX)) || (AngRadius1 < (N-pCentreY)))
                AngRadius = round(max([pCentreX, pCentreY, (M-pCentreX), (N-pCentreY)]));
            else
                AngRadius = AngRadius1;
            end
            polarheight= round((AngRadius-pRadius) * scale);
            x_iris = pCentreX;
            y_iris = pCentreY;
            [polar_iris] = convert_to_polar(thisFrame, x_iris, y_iris, AngRadius, pCentreX, pCentreY, pRadius,  polarheight, 360);
            
            % Shift IRIS
            polar_iris = shiftiris(polar_iris);
            
            % Localaize Iris Boundary
            [iboundary,  ~]= iris_boundary(polar_iris,pRadius);
            
            % Store this initial estimate of iris width
            irisInitialWidth(ii) = round(iboundary / scale);
            
        end % check defined pupil fit
    catch ME
        warning('Error processing frame %d',ii);
        disp(ME.message)
    end % try catch
end % loop through gray frames

% report completion of analysis
if strcmp(p.Results.verbosity,'full')
    toc
    fprintf('\n');
end

% Calculate the mean irisWidth across the frames
irisWidth = nanmean(irisInitialWidth);

% initialize variables to hold the results in the next parfor loop
irisFitData_X = nan(nFrames,1);
irisFitData_Y = nan(nFrames,1);
irisFitData_radius = nan(nFrames,1);
irisFitData_mask = zeros(videoSizeY,videoSizeX,nFrames,'uint8');

%% Now refine iris size and obtain iris mask

% alert the user
if strcmp(p.Results.verbosity,'full')
    tic
    fprintf(['Refining iris width and extracting mask. Started ' char(datetime('now')) '\n']);
    fprintf('| 0                      50                   100%% |\n');
    fprintf('.\n');
end

% loop over frames
parfor (ii = 1:nFrames, nWorkers)
    
    % increment the progress bar
    if strcmp(p.Results.verbosity,'full') && mod(ii,round(nFrames/50))==0
        fprintf('\b.\n');
    end
    
    % diagnostic try/catch
    try
        
        % check if there is a defined pupil fit. If so, proceed
        if ~isnan(pupilFitParams(ii,1))
            
            % get the frame
            thisFrame = squeeze(grayVideo(:,:,ii));
            
            % Calculate the pupil height from the corrected pupil perimeter
            pupil_height = min(find(~max(squeeze(pupilPerimeter.data(:,:,ii))')==0));
            
            % To maintain transparency with regard to the IrisSeg toolbox, we
            % transfer values to variable names used by those routines
            scale= 1;
            pCentreX=pupilFitParams(ii,1);
            pCentreY=pupilFitParams(ii,2);
            pRadius=sqrt(pupilFitParams(ii,3)/pi);
            
            %% Code pulled from irisseg_main.m, part of the IrisSeg toolbox
            % https://github.com/cdac-cvml/IrisSeg
            
            % Eyelid Occlusion Detection Module
            [eyelidmask, adaptImage,~, ~] = geteyelid(thisFrame, pCentreX, pCentreY, pRadius, irisWidth, pupil_height, scale);
            
            % Turn off warnings for nargchk being deprecated; we have
            % communicated to the IrisSeg folks that this needs to be fixed
            warningState = warning;
            warning('off','MATLAB:nargchk:deprecated');
            
            % Iris Boundary Refinement Module
            [final_CX, final_CY, Final_iRadius, irismask, ~] = iris_boundary_actual_double(thisFrame,adaptImage, eyelidmask, pCentreX, pCentreY, pRadius, irisWidth, scale);
            
            % Restore the warning state
            warning(warningState);
            
            % final adjustment for scale and store the results
            x_iris = final_CX / scale;
            y_iris = final_CY /scale;
            iRadius = Final_iRadius / scale;
            
            irisFitData_X(ii) = x_iris;
            irisFitData_Y(ii) = y_iris;
            irisFitData_radius(ii) = iRadius;
            irisFitData_mask(:,:,ii) = irismask;
            
        end % check defined pupil fit
    catch ME
        warning('Error processing frame %d',ii);
        disp(ME.message)
    end % try catch
end % loop through gray frames


%% Clean up and close

% gather the loop vars into the irisFitData structure
irisFitData.X = irisFitData_X;
irisFitData.Y = irisFitData_Y;
irisFitData.radius = irisFitData_radius;
irisFitData.mask=irisFitData_mask;

% save irisFitData
irisFitData.meta = p.Results;
save(irisFitFileName,'irisFitData');

% report completion of analysis
if strcmp(p.Results.verbosity,'full')
    toc
    fprintf('\n');
end

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

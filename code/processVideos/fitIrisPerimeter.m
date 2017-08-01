function fitIrisPerimeter(grayVideoName, perimeterFileName, pupilFileName, irisFileName, varargin)
% function fitIrisPerimeter(grayVideoName, perimeterFileName, pupilFileName, irisFileName, varargin)
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
%   perimeterFileName
%	pupilFileName - full path to the .mat file that contains fits
%       to the border of the pupil
%   irisFileName -
%
% Options (verbosity and display)
%   verbosity - controls console status updates
%
% Options (analysis)
% 	gammaCorrection - gamma correction to be applied to the video frames
%       (default 1, typical range [0.5 1.8])
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
p.addRequired('pupilFileName',@isstr);
p.addRequired('irisFileName',@isstr);

% Optional display params
p.addParameter('verbosity','none',@ischar);

% Optional analysis params
p.addParameter('irisGammaCorrection', .5, @isnumeric);
p.addParameter('irisTransparentLB',[0, 0, 10000, 0, -0.5*pi],@isnumeric);
p.addParameter('irisTransparentUB',[240,320,30000,0.3, 0.5*pi],@isnumeric);

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

%% Parse and check the parameters
p.parse(grayVideoName, perimeterFileName, pupilFileName, irisFileName, varargin{:})

nEllipseParams=5; % 5 params in the transparent ellipse form

if length(p.Results.irisTransparentLB)~=nEllipseParams
    error('Wrong number of elements in irisTransparentLB');
end
if length(p.Results.irisTransparentUB)~=nEllipseParams
    error('Wrong number of elements in irisTransparentUB');
end
if sum(p.Results.irisTransparentUB>=p.Results.irisTransparentLB)~=nEllipseParams
    error('Lower bounds must be equal to or less than upper bounds');
end


%% Prepare an anonymous function
% Create a non-linear constraint for the ellipse fit. If no parameters are
% given, then create an empty function handle (and thus have no non-linear
% constraint)
nonlinconst = @(x) restrictEccenByTheta(x,[p.Results.irisTransparentUB(4) p.Results.irisTransparentUB(4)]);


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
    thisFrame = imadjust(thisFrame,[],[],p.Results.irisGammaCorrection);
    grayVideo(:,:,ii) = rgb2gray (thisFrame);
end
% close the video object
clear videoInObj

% read in the pupil ellipse fit data
dataLoad = load(pupilFileName);
pupilData = dataLoad.pupilData;
clear dataLoad
pupilFitParams = pupilData.pPosteriorMeanTransparent;

% Load the pupil perimeter data. It will be a structure variable
% "perimeter", with the fields .data and .meta
dataLoad=load(perimeterFileName);
perimeter=dataLoad.perimeter;
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


%% Perform two passes on fit to iris perimeter

for pass=1:2
    
    % alert the user
    if strcmp(p.Results.verbosity,'full')
        tic
        switch pass
            case 1
                fprintf(['First pass iris fit. Started ' char(datetime('now')) '\n']);
            case 2
                fprintf(['Second pass iris fit. Started ' char(datetime('now')) '\n']);
        end
        fprintf('| 0                      50                   100%% |\n');
        fprintf('.\n');
    end
    
    % Create an anonymous function for ellipse fitting, varying the bounds
    % depending upon if we are on the first or second pass
    switch pass
        case 1
            lb = p.Results.irisTransparentLB;
            ub = p.Results.irisTransparentUB;
            obtainEllipseLikelihood = @(x,y) constrainedEllipseFit(x, y, ...
                lb, ub, nonlinconst);
        case 2
            lb = p.Results.irisTransparentLB;
            ub = p.Results.irisTransparentUB;
            medianIrisArea = nanmedian(irisData_pEllipseFitTransparent(:,3));
            lb(3) = medianIrisArea * 0.95;
            ub(3) = medianIrisArea * 1.05;
            obtainEllipseLikelihood = @(x,y) constrainedEllipseFit(x, y, ...
                lb, ub, nonlinconst);
    end
    
    % initialize variables to hold the results
    irisData_pEllipseFitTransparent = nan(nFrames,nEllipseParams);
    irisData_pEllipseFitHessianSD = nan(nFrames,nEllipseParams);
    
    
    % Loop through gray frames to find the iris
    parfor (ii = 1:nFrames, nWorkers)
        
        % increment the progress bar
        if strcmp(p.Results.verbosity,'full') && mod(ii,round(nFrames/50))==0
            fprintf('\b.\n');
        end
                    
            % check if there is a defined pupil fit. If so, proceed
            if ~isnan(pupilFitParams(ii,1))
                
                % get the frame
                thisFrame = squeeze(grayVideo(:,:,ii));
                
                % Calculate the pupil height from the corrected pupil perimeter
                pupil_height = find(~max(squeeze(perimeter.data(:,:,ii))')==0, 1 );
                
                % To maintain transparency with regard to the IrisSeg toolbox, we
                % transfer values to variable names used by those routines
                scale= 1;
                pCentreX=pupilFitParams(ii,1);
                pCentreY=pupilFitParams(ii,2);
                pRadius=sqrt(pupilFitParams(ii,3)/pi);
                N = videoSizeY;
                M = videoSizeX;
                
                try
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Code block from irisseg_main.m, part of the IrisSeg toolbox
                %   https://github.com/cdac-cvml/IrisSeg
                %
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
                
                % Obtain the iris width
                irisWidth = round(iboundary / scale);
                
                % Eyelid Occlusion Detection Module
                [eyelidmask, adaptImage,~, ~] = geteyelid(thisFrame, pCentreX, pCentreY, pRadius, irisWidth, pupil_height, scale);
                
                % Turn off warnings for nargchk being deprecated; we have
                % communicated to the IrisSeg folks that this needs to be fixed
                warningState = warning;
                warning('off','MATLAB:nargchk:deprecated');
                
                % Iris Boundary Refinement Module
                [~, ~, ~, irismask, ~] = iris_boundary_actual_double(thisFrame,adaptImage, eyelidmask, pCentreX, pCentreY, pRadius, irisWidth, scale);
                
                % Restore the warning state
                warning(warningState);
                %
                % END code block from IrisSeg
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                % Obtain the boundary of the irisMask
                bpTmp = bwboundaries(irismask);
                irisBoundary=bpTmp{1};
                nBoundaryPoints = size(irisBoundary,1);
                
                % Identify the sides of the boundary using the
                % slope of change in vertical position around the boundary
                movingSlopeSupportLength = 9;
                verticalChangeBoundaryIdx=find(abs(movingslope(irisBoundary(:,1), movingSlopeSupportLength)) > 0.75);
                if ~isempty(verticalChangeBoundaryIdx)
                    irisBoundary=irisBoundary(verticalChangeBoundaryIdx,:);
                    circleColor='red';
                end
                
                % fit an ellipse to the iris boundary points
                [pEllipseFitTransparent, pInitialFitHessianSD, ~] = ...
                    feval(obtainEllipseLikelihood,irisBoundary(:,2),irisBoundary(:,1));
                
                irisData_pEllipseFitTransparent(ii,:) = pEllipseFitTransparent';
                irisData_pEllipseFitHessianSD(ii,:) = pInitialFitHessianSD';
                catch
                    warning ('Error fitting the iris for frame %d', ii)
                    continue
                end
                
            end % check defined pupil fit
    end % loop through gray frames
    
    % report completion of this pass analysis
    if strcmp(p.Results.verbosity,'full')
        toc
        fprintf('\n');
    end
    
end % loop through first and second pass


%% Clean up and save
% gather the loop vars into the irisFitData structure
irisData.pEllipseFitTransparent=irisData_pEllipseFitTransparent;
irisData.pEllipseFitHessianSD=irisData_pEllipseFitHessianSD;

% save irisFitData
irisData.meta = p.Results;
save(irisFileName,'irisData');


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




function [c, ceq]=restrictEccenByTheta(transparentEllipseParams, constrainEccen_x_Theta)
% function [c, ceq]=restrictEccenByTheta(transparentEllipseParams,constrainEccen_x_Theta)
%
% This function implements a non-linear constraint upon the ellipse fit
% to the pupil boundary. The goal of the limit is to constrain theta to the
% cardinal axes, and more severely constrain eccentricity in the horizontal
% as compared to the vertical direction.

% First constraint (equality)
%  - the theta is on a cardinal axis (i.e., theta is from the set [-pi/2,0,pi/2])
ceq = mod(abs(transparentEllipseParams(5)),(pi/2));

% Second constraint (inequality)
%  - require the eccen to be less than constrainEccen_x_Theta, where this
%    has one value for horizontal ellipses (i.e., abs(theta)<pi/2) and a
%    second value for vertical ellipses.
if abs(transparentEllipseParams(5))<(pi/4)
    c = transparentEllipseParams(4) - constrainEccen_x_Theta(1);
else
    c = transparentEllipseParams(4) - constrainEccen_x_Theta(2);
end

end % restrictEccenByTheta func

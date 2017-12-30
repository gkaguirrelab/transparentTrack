function [pupilData] = fitPupilPerimeter(perimeterFileName, pupilFileName, varargin)
% Perform non-linear, constrained ellipse fitting to pupil perimeters
%
% Description:
%   This routine fits an ellipse to each frame of a video that contains the
%   perimeter of the pupil. A non-linear search routine is used, with
%   optional, non-linear constraints on the eccentricity and theta of the
%   ellipse that are informed by scene geometry. An estimate of the
%   standard deviation of the parameters of the best fitting ellipse are
%   calculated and stored as well.
%
% Notes:
%   Image coordinates -  MATLAB uses an "intrinsic" coordinate system such
%   that the center of each pixel in an image corresponds to its integer
%   indexed position. Thus a 3x3 pixel image in intrisic coordinates is
%   represented on a grid with xlim = [0.5 3.5] and ylim = [0.5 3.5], with
%   the origin being the top left corner of the image. This is done to
%   facilitate the handling of images in many of the built-in image
%   processing functions. This routine outputs the results in intrinsic
%   coordinates. Additinal information regarding the MATLAB image
%   coordinate system may be found here:
%       https://blogs.mathworks.com/steve/2013/08/28/introduction-to-spatial-referencing/
%
%   Parallel pool - Controlled by the key/value pair 'useParallel'. The
%   routine should gracefully fall-back on serial processing if the
%   parallel pool is unavailable. Each worker requires ~8 GB of memory to
%   operate. It is important to keep total RAM usage below the physical
%   memory limit to prevent swapping and a dramatic slow down in
%   processing. To use the parallel pool with TbTb, provide the identity of
%   the repo name in the 'tbtbRepoName', which is then used to configure
%   the workers.
%
% Inputs:
%   perimeterFileName     - Full path to a .mat file that contains the
%                           pupil perimeter data.
%   pupilFileName         - Full path to the .mat file in which to save
%                           the results of the ellipse fitting.
%
% Optional key/value pairs (display and I/O):
%  'verbosity'            - Level of verbosity. [none, full]
%
% Optional key/value pairs (flow control)
%  'nFrames'              - Analyze fewer than the total number of frames.
%  'useParallel'          - If set to true, use the Matlab parallel pool
%  'nWorkers'             - Specify the number of workers in the parallel
%                           pool. If undefined the default number will be
%                           used.
%  'tbtbProjectName'      - The workers in the parallel pool are configured
%                           by issuing a tbUseProject command for the
%                           project specified here.
%
% Optional key/value pairs (environment)
%  'tbSnapshot'           - This should contain the output of the
%                           tbDeploymentSnapshot performed upon the result
%                           of the tbUse command. This documents the state
%                           of the system at the time of analysis.
%  'timestamp'            - AUTOMATIC; The current time and date
%  'username'             - AUTOMATIC; The user
%  'hostname'             - AUTOMATIC; The host
%
% Optional key/value pairs (fitting)
%  'ellipseTransparentLB/UB' - Define the hard upper and lower boundaries
%                           for the ellipse fit, in units of pixels of the
%                           video. The default values selected here
%                           represent the physical and mathematical limits,
%                           as the constraint for the fit will be provided
%                           by the scene geometry. A mild constraint (0.75)
%                           is placed upon the eccentricity, corresponding
%                           to an aspect ration of 3:2.
%  'nSplits'              - The number of tests upon the spatial split-
%                           halves of the pupil boundary values to examine
%                           to estimate the SD of the fitting parameters.
%  'sceneGeometryFileName' - Full path to a sceneGeometry file. When the
%                           sceneGeometry is available, fitting is
%                           performed in terms of eye parameters instead of
%                           ellipse parameters
%
% Outputs:
%	pupilData             - A structure with multiple fields corresponding
%                           to the parameters, SDs, and errors of the
%                           fit. Different field names are used
%                           depending upon if a sceneGeometry constraint
%                           was or was not used.
%

%% Parse vargin for options passed here
p = inputParser; p.KeepUnmatched = true;

% Required
p.addRequired('perimeterFileName',@ischar);
p.addRequired('pupilFileName',@ischar);

% Optional display and I/O params
p.addParameter('verbosity','none',@ischar);

% Optional flow control params
p.addParameter('nFrames',Inf,@isnumeric);
p.addParameter('useParallel',false,@islogical);
p.addParameter('nWorkers',[],@(x)(isempty(x) | isnumeric(x)));
p.addParameter('tbtbRepoName','transparentTrack',@ischar);

% Optional environment params
p.addParameter('tbSnapshot',[],@(x)(isempty(x) | isstruct(x)));
p.addParameter('timestamp',char(datetime('now')),@ischar);
p.addParameter('hostname',char(java.lang.System.getProperty('user.name')),@ischar);
p.addParameter('username',char(java.net.InetAddress.getLocalHost.getHostName),@ischar);

% Optional analysis params
p.addParameter('ellipseTransparentLB',[0, 0, 800, 0, 0],@(x)(isempty(x) | isnumeric(x)));
p.addParameter('ellipseTransparentUB',[640,480,20000,0.75, pi],@(x)(isempty(x) | isnumeric(x)));
p.addParameter('nSplits',8,@isnumeric);
p.addParameter('sceneGeometryFileName',[],@(x)(isempty(x) | ischar(x)));
p.addParameter('ellipseFitLabel',[],@(x)(isempty(x) | ischar(x)));


%% Parse and check the parameters
p.parse(perimeterFileName, pupilFileName, varargin{:});

nEllipseParams=5; % 5 params in the transparent ellipse form
nEyeParams=3; % 3 values (azimuth, elevation, pupil radius) for eyeParams

%% Prepare some anonymous functions
% Create an anonymous function to return a rotation matrix given theta in
% radians
returnRotMat = @(theta) [cos(theta) -sin(theta); sin(theta) cos(theta)];

%% Load the pupil perimeter data and optionally sceneGeometry
% Load the pupil perimeter data. It will be a structure variable
% "perimeter", with the fields .data and .meta
dataLoad=load(perimeterFileName);
perimeter=dataLoad.perimeter;
clear dataLoad

if isempty(p.Results.sceneGeometryFileName)
    sceneGeometry=[];
else
    % load the sceneGeometry structure
    dataLoad=load(p.Results.sceneGeometryFileName);
    sceneGeometry=dataLoad.sceneGeometry;
    clear dataLoad
end

% determine how many frames we will process
if p.Results.nFrames == Inf
    nFrames=size(perimeter.data,1);
else
    nFrames = p.Results.nFrames;
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


%% Calculate an ellipse fit for each video frame

% Recast perimeter into a sliced cell array to reduce par for
% broadcast overhead
frameCellArray = perimeter.data(1:nFrames);
clear perimeter

% Set-up other variables to be non-broadcast
verbosity = p.Results.verbosity;
ellipseTransparentLB = p.Results.ellipseTransparentLB;
ellipseTransparentUB = p.Results.ellipseTransparentUB;
nSplits = p.Results.nSplits;

% Alert the user
if strcmp(p.Results.verbosity,'full')
    tic
    fprintf(['Ellipse fitting to pupil perimeter. Started ' char(datetime('now')) '\n']);
    fprintf('| 0                      50                   100%% |\n');
    fprintf('.\n');
end

% Loop through the frames
parfor (ii = 1:nFrames, nWorkers)
    
    % Update progress
    if strcmp(verbosity,'full')
        if mod(ii,round(nFrames/50))==0
            fprintf('\b.\n');
        end
    end
    
    % Initialize the results variables
    ellipseParamsTransparent=NaN(1,nEllipseParams);
    ellipseParamsSplitsSD=NaN(1,nEllipseParams);
    ellipseParamsObjectiveError=NaN(1);
    eyeParams=NaN(1,nEyeParams);
    eyeParamsSplitsSD=NaN(1,nEyeParams);
    eyeParamsObjectiveError=NaN(1);
    pFitTransparentSplit=NaN(1,nSplits,nEllipseParams);
    pFitEyeParamSplit=NaN(1,nSplits,nEyeParams);
    
    % get the boundary points
    Xp = frameCellArray{ii}.Xp;
    Yp = frameCellArray{ii}.Yp;
    
    % fit an ellipse to the boundary (if any points exist)
    if ~isempty(Xp) && ~isempty(Yp)
        try % this is to have information on which frame caused an error
            
            % Obtain the fit to the veridical data
            if isempty(sceneGeometry)
                [ellipseParamsTransparent, ellipseParamsObjectiveError] = ...
                    constrainedEllipseFit(Xp, Yp, ...
                    ellipseTransparentLB, ...
                    ellipseTransparentUB, ...
                    []);
            else
                [eyeParams, eyeParamsObjectiveError] = eyeParamEllipseFit(Xp, Yp, sceneGeometry);
                ellipseParamsTransparent = pupilProjection_fwd(eyeParams, sceneGeometry);
            end
            
            % Re-calculate fit for splits of data points, if requested
            if nSplits == 0
                if isempty(sceneGeometry)
                    ellipseParamsSplitsSD=NaN(1,nEllipseParams);
                else
                    eyeParamsSplitsSD=NaN(1,nEyeParams);
                end
            else
                % Find the center of the pupil boundary points, place the boundary
                % points in a matrix and shift them to the center position
                xCenter=mean(Xp); yCenter=mean(Yp);
                centerMatrix = repmat([xCenter'; yCenter'], 1, length(Xp));
                
                % Prepare a variable to hold the results of the split data
                % fits
                if isempty(sceneGeometry)
                    pFitTransparentSplit=NaN(2,nSplits,nEllipseParams);
                else
                    pFitEyeParamSplit=NaN(2,nSplits,nEyeParams);
                end
                % Loop across the number of requested splits
                for ss=1:nSplits
                    % Rotate the data and split in half through the center
                    theta=((pi/2)/nSplits)*ss;
                    forwardPoints = feval(returnRotMat,theta) * ([Xp,Yp]' - centerMatrix) + centerMatrix;
                    splitIdx1 = find((forwardPoints(1,:) < median(forwardPoints(1,:))))';
                    splitIdx2 = find((forwardPoints(1,:) >= median(forwardPoints(1,:))))';
                    % Fit the split sets of pupil boundary points
                    if isempty(sceneGeometry)
                        pFitTransparentSplit(1,ss,:) = ...
                            constrainedEllipseFit(Xp(splitIdx1), Yp(splitIdx1), ...
                            ellipseTransparentLB, ...
                            ellipseTransparentUB, ...
                            []);
                        pFitTransparentSplit(2,ss,:) = ...
                            constrainedEllipseFit(Xp(splitIdx2), Yp(splitIdx2), ...
                            ellipseTransparentLB, ...
                            ellipseTransparentUB, ...
                            []);
                    else
                        pFitEyeParamSplit(1,ss,:) = eyeParamEllipseFit(Xp(splitIdx1), Yp(splitIdx1), sceneGeometry);
                        pFitEyeParamSplit(2,ss,:) = eyeParamEllipseFit(Xp(splitIdx1), Yp(splitIdx1), sceneGeometry);
                    end
                end % loop through splits
                
                % Calculate the SD of the parameters across splits
                if isempty(sceneGeometry)
                    ellipseParamsSplitsSD=nanstd(reshape(pFitTransparentSplit,ss*2,nEllipseParams));
                else
                    eyeParamsSplitsSD=nanstd(reshape(pFitEyeParamSplit,ss*2,nEyeParams));
                end
            end % check if we want to do splits
            
        catch ME
            warning ('Error while processing frame: %d', ii)
        end % try catch
    end % check if there are pupil boundary data to be fit
    
    % store results
    loopVar_ellipseParamsTransparent(ii,:) = ellipseParamsTransparent';
    loopVar_ellipseParamsSplitsSD(ii,:) = ellipseParamsSplitsSD';
    loopVar_ellipseParamsObjectiveError(ii) = ellipseParamsObjectiveError;
    if ~isempty(sceneGeometry)
        loopVar_eyeParams(ii,:) = eyeParams';
        loopVar_eyeParamsSplitsSD(ii,:) = eyeParamsSplitsSD';
        loopVar_eyeParamsObjectiveError(ii) = eyeParamsObjectiveError;
    end
end % loop over frames

% alert the user that we are done with the fit loop
if strcmp(p.Results.verbosity,'full')
    toc
    fprintf('\n');
end

%% Clean up and save

% Check if the pupilData file already exists. If so, load it. This allows
% us to preserve the Unconstrained results when a subsequent, sceneGeometry
% constrained fit is constructed.
if exist(p.Results.pupilFileName, 'file') == 2
    dataLoad=load(pupilFileName);
    pupilData=dataLoad.pupilData;
    clear dataLoad
end

% Establish a label to save the fields of the ellipse fit data
if isempty(p.Results.ellipseFitLabel)
    if isempty(sceneGeometry)
        ellipseFitLabel = 'initial';
    else
        ellipseFitLabel = 'sceneConstrained';
    end
else
    ellipseFitLabel = p.Results.ellipseFitLabel;
end

% Store the ellipse fit data in informative fields and add meta data
pupilData.(ellipseFitLabel).ellipse.values = loopVar_ellipseParamsTransparent;
if isempty(sceneGeometry)
    pupilData.(ellipseFitLabel).ellipse.RMSE = loopVar_ellipseParamsObjectiveError';
else
    pupilData.(ellipseFitLabel).ellipse.RMSE = loopVar_eyeParamsObjectiveError';
end
if nSplits~=0 && isempty(sceneGeometry)
    pupilData.(ellipseFitLabel).ellipse.splitsSD = loopVar_ellipseParamsSplitsSD;
end
pupilData.(ellipseFitLabel).ellipse.meta.ellipseForm = 'transparent';
pupilData.(ellipseFitLabel).ellipse.meta.labels = {'x','y','area','eccentricity','theta'};
pupilData.(ellipseFitLabel).ellipse.meta.units = {'pixels','pixels','squared pixels','non-linear eccentricity','rads'};
pupilData.(ellipseFitLabel).ellipse.meta.coordinateSystem = 'intrinsic image';
if ~isempty(sceneGeometry)
    pupilData.(ellipseFitLabel).eyeParams.values = loopVar_eyeParams;
    if nSplits~=0
        pupilData.(ellipseFitLabel).eyeParams.splitsSD = loopVar_eyeParamsSplitsSD;
    end
    pupilData.(ellipseFitLabel).eyeParams.meta.labels = {'azimuth','elevation','pupil radius'};
    pupilData.(ellipseFitLabel).eyeParams.meta.units = {'deg','deg','mm'};
    pupilData.(ellipseFitLabel).eyeParams.meta.coordinateSystem = 'head fixed (extrinsic)';
end
pupilData.(ellipseFitLabel).meta = p.Results;
% save the ellipse fit results
save(p.Results.pupilFileName,'pupilData')


%% Delete the parallel pool
if p.Results.useParallel
    if strcmp(p.Results.verbosity,'full')
        tic
        fprintf(['Closing parallel pool. Started ' char(datetime('now')) '\n']);
    end
    poolObj = gcp;
    if ~isempty(poolObj)
        delete(poolObj);
    end
    if strcmp(p.Results.verbosity,'full')
        toc
        fprintf('\n');
    end
end


end % function


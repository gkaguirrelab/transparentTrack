function [pupilData] = fitPupilPerimeter(perimeterFileName, pupilFileName, varargin)
% Perform non-linear, constrained ellipse fitting to pupil perimeters
%
% Syntax:
%  [pupilData] = fitPupilPerimeter(perimeterFileName, pupilFileName)
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
%   processing functions. This routine outputs results in intrinsic
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
%                           by the scene geometry. A mild constraint (0.6)
%                           is placed upon the eccentricity, corresponding
%                           to an aspect ration of 4:5.
%  'eyePoseLB/UB'         - Upper and lower bounds on the eyePose
%                           [azimuth, elevation, pupil radius]. Biological
%                           limits in eye rotation and pupil size would
%                           suggest boundaries of [±35, ±25, 0.5-5]. Note,
%                           however, that these angles are relative to the
%                           center of projection, not the primary position
%                           of the eye. Therefore, in circumstances in
%                           which the camera is viewing the eye from an
%                           off-center angle, the bounds will need to be
%                           shifted accordingly.
%  'nSplits'              - The number of tests upon the spatial split-
%                           halves of the pupil boundary values to examine
%                           to estimate the SD of the fitting parameters.
%  'sceneGeometryFileName' - Full path to a sceneGeometry file. When the
%                           sceneGeometry is available, fitting is
%                           performed in terms of eye parameters instead of
%                           ellipse parameters
%  'fitLabel'             - The field name in the pupilData structure where
%                           the results of the fitting will be stored.
%  'useRayTracing'        - Logical; default false. Using ray tracing in
%                           the camera translation search improves accuracy
%                           slightly, but increases search time by about
%                           25x.

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
p.addParameter('ellipseTransparentLB',[0,0,800,0,0],@(x)(isempty(x) | isnumeric(x)));
p.addParameter('ellipseTransparentUB',[640,480,20000,0.6,pi],@(x)(isempty(x) | isnumeric(x)));
p.addParameter('eyePoseLB',[-35,-25,0,0.25],@(x)(isempty(x) | isnumeric(x)));
p.addParameter('eyePoseUB',[35,25,0,4],@(x)(isempty(x) | isnumeric(x)));
p.addParameter('nSplits',2,@isnumeric);
p.addParameter('sceneGeometryFileName',[],@(x)(isempty(x) | ischar(x)));
p.addParameter('fitLabel',[],@(x)(isempty(x) | ischar(x)));
p.addParameter('useRayTracing',false,@islogical);


%% Parse and check the parameters
p.parse(perimeterFileName, pupilFileName, varargin{:});

nEllipseParams=5; % 5 params in the transparent ellipse form
nEyePoseParams=4; % 4 values (azimuth, elevation, torsion, pupil radius) for eyePoses


%% Load data
% Load the pupil perimeter data. It will be a structure variable
% "perimeter", with the fields .data and .meta
dataLoad=load(perimeterFileName);
perimeter=dataLoad.perimeter;
clear dataLoad

% Optionally load a sceneGeometry file
if isempty(p.Results.sceneGeometryFileName)
    sceneGeometry=[];
else
    % load the sceneGeometry structure
    dataLoad=load(p.Results.sceneGeometryFileName);
    sceneGeometry=dataLoad.sceneGeometry;
    clear dataLoad
end

% If sceneGeometry is defined, prepare the ray tracing functions
if ~isempty(sceneGeometry) && p.Results.useRayTracing
    if strcmp(p.Results.verbosity,'full')
        fprintf('Assembling ray tracing functions.\n');
    end
    [rayTraceFuncs] = assembleRayTraceFuncs( sceneGeometry );
else
    rayTraceFuncs = [];
end

% Optionally load the pupilData file
if exist(p.Results.pupilFileName, 'file') == 2
    dataLoad=load(pupilFileName);
    pupilData=dataLoad.pupilData;
    clear dataLoad
else
    pupilData=[];
end

% determine how many frames we will process
if p.Results.nFrames == Inf
    nFrames=size(perimeter.data,1);
else
    nFrames = p.Results.nFrames;
end


%% Prepare some functions
% Create an anonymous function to return a rotation matrix given theta in
% radians
returnRotMat = @(theta) [cos(theta) -sin(theta); sin(theta) cos(theta)];


%% Set up the parallel pool
if p.Results.useParallel
    nWorkers = startParpool( p.Results.nWorkers, p.Results.tbtbRepoName, p.Results.verbosity );
else
    nWorkers=0;
end


%% Calculate an ellipse fit for each video frame

% Recast perimeter into a sliced cell array to reduce parfor
% broadcast overhead
frameCellArray = perimeter.data(1:nFrames);
clear perimeter

% Set-up other variables to be non-broadcast
verbosity = p.Results.verbosity;
ellipseTransparentLB = p.Results.ellipseTransparentLB;
ellipseTransparentUB = p.Results.ellipseTransparentUB;
eyePoseLB = p.Results.eyePoseLB;
eyePoseUB = p.Results.eyePoseUB;
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
%for ii = 311:nFrames
    
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
    eyePose=NaN(1,nEyePoseParams);
    eyePoseSplitsSD=NaN(1,nEyePoseParams);
    eyePoseObjectiveError=NaN(1);
    pFitTransparentSplit=NaN(1,nSplits,nEllipseParams);
    pFitEyePoseSplit=NaN(1,nSplits,nEyePoseParams);
    
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
                % Identify the best fitting eye parameters for the  the
                % pupil perimeter
                [eyePose, eyePoseObjectiveError] = ...
                    eyePoseEllipseFit(Xp, Yp, sceneGeometry, rayTraceFuncs, 'eyePoseLB', eyePoseLB, 'eyePoseUB', eyePoseUB);
                % Obtain the parameters of the ellipse
                ellipseParamsTransparent = ...
                    pupilProjection_fwd(eyePose, sceneGeometry, rayTraceFuncs);
            end
            
            % Re-calculate fit for splits of data points, if requested
            if nSplits == 0
                if isempty(sceneGeometry)
                    ellipseParamsSplitsSD=NaN(1,nEllipseParams);
                else
                    eyePoseSplitsSD=NaN(1,nEyePoseParams);
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
                    pFitEyePoseSplit=NaN(2,nSplits,nEyePoseParams);
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
                        % We don't have sceneGeometry defined, so fit an
                        % ellipse to the splits of the pupil perimeter
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
                        % We do have sceneGeometry, so search for eyePose
                        % that best fit the splits of the pupil perimeter.
                        % To speed up the search, we do not use ray tracing
                        % here, as we are not interested in the absolute
                        % values of the fits, but instead just their
                        % variation.
                        pFitEyePoseSplit(1,ss,:) = ...
                            eyePoseEllipseFit(Xp(splitIdx1), Yp(splitIdx1), sceneGeometry, [], 'x0', eyePose, 'eyePoseLB', eyePoseLB, 'eyePoseUB', eyePoseUB);
                        pFitEyePoseSplit(2,ss,:) = ...
                            eyePoseEllipseFit(Xp(splitIdx2), Yp(splitIdx2), sceneGeometry, [], 'x0', eyePose, 'eyePoseLB', eyePoseLB, 'eyePoseUB', eyePoseUB);
                        % Obtain the ellipse parameeters that correspond
                        % the eyePose
                        pFitTransparentSplit(1,ss,:) = ...
                            pupilProjection_fwd(pFitEyePoseSplit(1,ss,:), sceneGeometry, []);
                        pFitTransparentSplit(2,ss,:) = ...
                            pupilProjection_fwd(pFitEyePoseSplit(2,ss,:), sceneGeometry, []);
                    end
                end % loop through splits
                
                % Calculate the SD of the parameters across splits
                ellipseParamsSplitsSD=nanstd(reshape(pFitTransparentSplit,ss*2,nEllipseParams));
                if ~isempty(sceneGeometry)
                    eyePoseSplitsSD=nanstd(reshape(pFitEyePoseSplit,ss*2,nEyePoseParams));
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
        loopVar_eyePoses(ii,:) = eyePose';
        loopVar_eyePosesSplitsSD(ii,:) = eyePoseSplitsSD';
        loopVar_eyePosesObjectiveError(ii) = eyePoseObjectiveError;
    end
    
end % loop over frames

% alert the user that we are done with the fit loop
if strcmp(p.Results.verbosity,'full')
    toc
    fprintf('\n');
end

%% Clean up and save

% Establish a label to save the fields of the ellipse fit data
if isempty(p.Results.fitLabel)
    if isempty(sceneGeometry)
        fitLabel = 'initial';
    else
        fitLabel = 'sceneConstrained';
    end
else
    fitLabel = p.Results.fitLabel;
end

% Store the ellipse fit data in informative fields
pupilData.(fitLabel).ellipses.values = loopVar_ellipseParamsTransparent;
if isempty(sceneGeometry)
    pupilData.(fitLabel).ellipses.RMSE = loopVar_ellipseParamsObjectiveError';
else
    pupilData.(fitLabel).ellipses.RMSE = loopVar_eyePosesObjectiveError';
end
if nSplits~=0
    pupilData.(fitLabel).ellipses.splitsSD = loopVar_ellipseParamsSplitsSD;
end
pupilData.(fitLabel).ellipses.meta.ellipseForm = 'transparent';
pupilData.(fitLabel).ellipses.meta.labels = {'x','y','area','eccentricity','theta'};
pupilData.(fitLabel).ellipses.meta.units = {'pixels','pixels','squared pixels','non-linear eccentricity','rads'};
pupilData.(fitLabel).ellipses.meta.coordinateSystem = 'intrinsic image';
if ~isempty(sceneGeometry)
    pupilData.(fitLabel).eyePoses.values = loopVar_eyePoses;
    if nSplits~=0
        pupilData.(fitLabel).eyePoses.splitsSD = loopVar_eyePosesSplitsSD;
    end
    pupilData.(fitLabel).eyePoses.meta.labels = {'azimuth','elevation','torsion','pupil radius'};
    pupilData.(fitLabel).eyePoses.meta.units = {'deg','deg','deg','mm'};
    pupilData.(fitLabel).eyePoses.meta.coordinateSystem = 'head fixed (extrinsic)';
end

% add meta data
pupilData.(fitLabel).meta = p.Results;

% save the ellipse fit results
save(p.Results.pupilFileName,'pupilData')

end % function


function [pupilData] = fitPupilPerimeter(perimeterFileName, pupilFileName, varargin)
% fitPupilPerimeter(perimeterVideoFileName, varargin)
%
% This routine fits an ellipse to each frame of a video that contains the
% perimeter of the pupil. An ellipse is fit to each frame of the video
% using a non-linear search routine, with some constraints on size and
% aspect ratio of the solution. An estimate of the standard deviation of
% the parameters of the best fitting ellipse is stored as well.
%
% A note on ellipse parameterization: an ellipse can be specified in
% multiple forms. Within the context of this routine, and in saved files,
% ellipses are considered in "transparent" form (our coinage):
%
%   center (cx,cy), area (a), eccentricity (e), angle of tilt (theta)
%
% We use this parameterization to allow us to constrain fits with regard to
% these values (specifically area and eccentricity).
%
% NOTES REGARDING image coordinates
%
% Matlab uses an "intrinsic" coordinate system such that the center of each
% pixel in an image corresponds to its integer indexed position. Thus a 3x3
% pixel image in intrisic coordinates is represented on a grid with xlim =
% [0.5 3.5] and ylim = [0.5 3.5], with the origin being the top left corner
% of the image. This is done to facilitate the handling of images in many
% of the built-in image processing functions. This routine outputs the
% results in intrinsic coordinates. Additinal information regarding the
% Matlab image coordinate system may be found here:
%   https://blogs.mathworks.com/steve/2013/08/28/introduction-to-spatial-referencing/
%
% NOTES REGARDING USE OF PARALLEL POOL
%
% The initial ellipse fitting is conducted within a parfor loop. The
% parallel pool will not be used unless the key/value pair 'useParallel' is
% set to true. The routine should gracefully fall-back on serial processing
% if the parallel pool is unavailable.
%
% Each worker requires ~8 GB of memory to operate. It is important to keep
% total RAM usage below the physical memory limit to prevent swapping and a
% dramatic slow down in processing.
%
% To use the parallel pool with TbTb, provide the identity of the repo name
% in the 'tbtbRepoName', which is then used to configure the workers.
%
% INPUTS:
%   perimeterFileName: full path to a .mat file that contains the perimeter
%     data varaible. Points on the boundary of the pupil should have a
%     value of unity, and the frame should be otherwise zero-filled. A
%     frame that has no information regarding the pupil (e.g., during a
%     blink) should be zero-filled.
%   pupilFileName: full path to the .mat file in which to save
%     pupil tracking information.
%
% Optional key/value pairs (verbosity)
%  'verbosity' - level of verbosity. [none, full]
%
% Optional key/value pairs (flow control)
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
% Optional key/value pairs (fitting parameters)
%
%  'ellipseTransparentLB/UB' - Define the hard upper and lower boundaries
%     for the ellipse fit, in units of pixels of the video. The center
%     points should be constrained to the size of the video.
%     Eccentricity is related to ratio of the semimajor and semiminor axes,
%     and can be calculated using:
%           eccentricity = axes2ecc(semimajor, semiminor)
%     For example, if we wish to prevent ellipses with an aspect ratio
%     greater than 3 : 2, this gives us an eccentricity UB of ~0.75.
%   'nSplits' - The number of tests upon the spatial split-halves of the
%     pupil boundary values to examine to estimate a likelihood SD.
%   'nBoots' - The number of bootstrap resamples of the pupil boundary
%     points to perform to estimate a likelihood SD. This was found in
%     testing to not be useful, as the pupil boundary is oversampled, so
%     this is left at a default of zero.
%
% OUTPUTS:
%   pupilData: A structure with multiple fields corresponding to the
%     parameters, SDs, and errors of the initial and final ellipse fits.


%% Parse vargin for options passed here
p = inputParser; p.KeepUnmatched = true;

% Required
p.addRequired('perimeterFileName',@ischar);
p.addRequired('pupilFileName',@ischar);

% Optional display and I/O params
p.addParameter('verbosity','none',@ischar);

% Optional fitting params
p.addParameter('ellipseTransparentLB',[0, 0, 800, 0, 0],@(x)(isempty(x) | isnumeric(x)));
p.addParameter('ellipseTransparentUB',[640,480,20000,0.75, pi],@(x)(isempty(x) | isnumeric(x)));
p.addParameter('nSplits',8,@isnumeric);

% Optional analysis params -- sceneGeometry fitting constraint
p.addParameter('sceneGeometryFileName',[],@(x)(isempty(x) | ischar(x)));

% Optional flow control params
p.addParameter('nFrames',Inf,@isnumeric);
p.addParameter('useParallel',false,@islogical);
p.addParameter('nWorkers',[],@(x)(isempty(x) | isnumeric(x)));
p.addParameter('tbtbRepoName','transparentTrack',@ischar);

% Environment parameters
p.addParameter('tbSnapshot',[],@(x)(isempty(x) | isstruct(x)));
p.addParameter('timestamp',char(datetime('now')),@ischar);
p.addParameter('hostname',char(java.lang.System.getProperty('user.name')),@ischar);
p.addParameter('username',char(java.net.InetAddress.getLocalHost.getHostName),@ischar);

%% Parse and check the parameters
p.parse(perimeterFileName, pupilFileName, varargin{:});

nEllipseParams=5; % 5 params in the transparent ellipse form

%% Prepare some anonymous functions and load the pupil perimeter data
% Create a non-linear constraint for the ellipse fit. If no parameters are
% given, then create an empty function handle (and thus have no non-linear
% constraint)
if isempty(p.Results.sceneGeometryFileName)
    nonlinconst = [];
else
    % load the sceneGeometry structure
    dataLoad=load(p.Results.sceneGeometryFileName);
    sceneGeometry=dataLoad.sceneGeometry;
    clear dataLoad
    
    nonlinconst = @(transparentEllipseParams) constrainEllipseBySceneGeometry(...
        transparentEllipseParams, ...
        sceneGeometry);
end

% Create an anonymous function to return a rotation matrix given theta in
% radians
returnRotMat = @(theta) [cos(theta) -sin(theta); sin(theta) cos(theta)];

% Load the pupil perimeter data. It will be a structure variable
% "perimeter", with the fields .data and .meta
dataLoad=load(perimeterFileName);
perimeter=dataLoad.perimeter;
clear dataLoad

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
    ellipseParamsConstraintError=NaN(1);
    
    % get the boundary points
    Xp = frameCellArray{ii}.Xp;
    Yp = frameCellArray{ii}.Yp;
    
    % fit an ellipse to the boundary (if any points exist)
    if ~isempty(Xp) && ~isempty(Yp)
        try % this is to have information on which frame caused an error
            
            % Obtain the fit to the veridical data
            [ellipseParamsTransparent, ellipseParamsObjectiveError, ellipseParamsConstraintError] = ...
                constrainedEllipseFit(Xp, Yp, ...
                ellipseTransparentLB, ...
                ellipseTransparentUB, ...
                nonlinconst);           
                        
            % Re-calculate fit for splits of data points, if requested
            if nSplits == 0
                ellipseParamsSplitsSD=NaN(1,nEllipseParams);
            else
                % Find the center of the pupil boundary points, place the boundary
                % points in a matrix and shift them to the center position
                xCenter=mean(Xp); yCenter=mean(Yp);
                centerMatrix = repmat([xCenter'; yCenter'], 1, length(Xp));
                
                % Rotate the data and split in half through the center
                pFitTransparentSplit=NaN(2,nSplits,nEllipseParams);
                for ss=1:nSplits
                    theta=((pi/2)/nSplits)*ss;
                    forwardPoints = feval(returnRotMat,theta) * ([Xp,Yp]' - centerMatrix) + centerMatrix;
                    splitIdx1 = find((forwardPoints(1,:) < median(forwardPoints(1,:))))';
                    splitIdx2 = find((forwardPoints(1,:) >= median(forwardPoints(1,:))))';
                    pFitTransparentSplit(1,ss,:) = ...
                        constrainedEllipseFit(Xp(splitIdx1), Yp(splitIdx1), ...
                        ellipseTransparentLB, ...
                        ellipseTransparentUB, ...
                        nonlinconst);
                    pFitTransparentSplit(2,ss,:) = ...
                        constrainedEllipseFit(Xp(splitIdx2), Yp(splitIdx2), ...
                        ellipseTransparentLB, ...
                        ellipseTransparentUB, ...
                        nonlinconst);
                end % loop through splits
                
                % Calculate the SD of the parameters across splits
                ellipseParamsSplitsSD=nanstd(reshape(pFitTransparentSplit,ss*2,nEllipseParams));
            end % check if we want to do splits
            
        catch ME
            warning ('Error while processing frame: %d', ii)
        end % try catch
    end % check if there are pupil boundary data to be fit
    
    % store results
    loopVar_ellipseParamsTransparent(ii,:) = ellipseParamsTransparent';
    loopVar_ellipseParamsSplitsSD(ii,:) = ellipseParamsSplitsSD';
    loopVar_ellipseParamsObjectiveError(ii) = ellipseParamsObjectiveError;
    loopVar_ellipseParamsConstraintError(ii) = ellipseParamsConstraintError;
    
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

% Store the ellipse fit data in informative fields and add meta data
if isempty(p.Results.sceneGeometryFileName)
    pupilData.meta.fitPupilPerimeterUnconstrained = p.Results;
    pupilData.meta.fitPupilPerimeterUnconstrained.coordinateSystem = 'intrinsicCoordinates(pixels)';
    pupilData.ellipseParamsUnconstrained_mean = loopVar_ellipseParamsTransparent;
    pupilData.ellipseParamsUnconstrained_rmse = loopVar_ellipseParamsObjectiveError';
    pupilData.ellipseParamsUnconstrained_constraint = loopVar_ellipseParamsConstraintError';
    if nSplits~=0
        pupilData.ellipseParamsUnconstrained_splitsSD = loopVar_ellipseParamsSplitsSD;
    end
else
    pupilData.meta.fitPupilPerimeterSceneConstrained = p.Results;
    pupilData.meta.fitPupilPerimeterSceneConstrained.coordinateSystem = 'intrinsicCoordinates(pixels)';
    pupilData.ellipseParamsSceneConstrained_mean = loopVar_ellipseParamsTransparent;
    pupilData.ellipseParamsSceneConstrained_rmse = loopVar_ellipseParamsObjectiveError';
    pupilData.ellipseParamsSceneConstrained_constraint = loopVar_ellipseParamsConstraintError';
    if nSplits~=0
        pupilData.ellipseParamsSceneConstrained_splitsSD = loopVar_ellipseParamsSplitsSD;
    end
end

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


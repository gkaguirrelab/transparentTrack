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
% Optional key/value pairs (fitting parameters)
%
%  'ellipseTransparentLB/UB' - Define the hard upper and lower boundaries
%     for the ellipse fit, in units of pixels of the video. The center
%     points should be constrained to the size of the video. Eccentricity
%     is related to ratio of the semimajor and semiminor axes, and can be
%     calculated using:
%           eccentricity = axes2ecc(semimajor, semiminor)
%     If we wish to prevent ellipses with an aspect ratio greater than
%     1.1 : 1, this gives us an eccentricity UB threshold of ~0.417.
%   'constrainEccen_x_Theta' - If defined, the ellipse fitting will be
%     constrained to allow only eccentric ellipses aligned with the
%     vertical or horizontal axes. Further, the two values provided will
%     differently limit the eccentricity of ellipses on the horizontal and
%     vertical axes, respectively.
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
p.addParameter('ellipseTransparentLB',[0, 0, 800, 0, -0.5*pi],@isnumeric);
p.addParameter('ellipseTransparentUB',[640,480,20000,0.5, 0.5*pi],@isnumeric);
p.addParameter('nSplits',8,@isnumeric);
p.addParameter('nBoots',0,@isnumeric);

% Optional analysis params -- sceneGeometry fitting constraint
p.addParameter('sceneGeometryFileName',[],@(x)(isempty(x) | ischar(x)));
p.addParameter('constraintMarginEccenMultiplier',1.05,@isnumeric);
p.addParameter('constraintMarginThetaDegrees',5,@isnumeric);

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

if length(p.Results.ellipseTransparentLB)~=nEllipseParams
    error('Wrong number of elements in ellipseTransparentLB');
end
if length(p.Results.ellipseTransparentUB)~=nEllipseParams
    error('Wrong number of elements in ellipseTransparentUB');
end
if sum(p.Results.ellipseTransparentUB>=p.Results.ellipseTransparentLB)~=nEllipseParams
    error('Lower bounds must be equal to or less than upper bounds');
end


%% Prepare some anonymous functions and load the pupil perimeter data
% Create a non-linear constraint for the ellipse fit. If no parameters are
% given, then create an empty function handle (and thus have no non-linear
% constraint)
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
        sceneGeometry, ...
        p.Results.constraintMarginEccenMultiplier, ...
        p.Results.constraintMarginThetaDegrees);
end

% Create an anonymous function for ellipse fitting
obtainPupilLikelihood = @(Xp,Yp) constrainedEllipseFit(Xp, Yp, ...
    p.Results.ellipseTransparentLB, p.Results.ellipseTransparentUB, nonlinconst);

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
    nFrames=size(perimeter.data,3);
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
    if strcmp(p.Results.verbosity,'full')
        if mod(ii,round(nFrames/50))==0
            fprintf('\b.\n');
        end
    end
    try % this is to have information on which frame caused an error
        % get the data frame
        thisFrame = squeeze(perimeter.data(:,:,ii));
        
        % get the boundary points
        [Yc, Xc] = ind2sub(size(thisFrame),find(thisFrame));
        
        % fit an ellipse to the boundary (if any points exist)
        if isempty(Xc) || isempty(Yc)
            pInitialFitTransparent=NaN(1,nEllipseParams);
            pInitialFitHessianSD=NaN(1,nEllipseParams);
            pInitialFitSplitsSD=NaN(1,nEllipseParams);
            pInitialFitBootsSD=NaN(1,nEllipseParams);
            pInitialFitError=NaN(1);
        else
            % Obtain the fit to the veridical data
            [pInitialFitTransparent, pInitialFitHessianSD, pInitialFitError] = ...
                feval(obtainPupilLikelihood,Xc, Yc);

            % Re-calculate fit for splits of data points, if requested
            if p.Results.nSplits == 0
                pInitialFitSplitsSD=NaN(1,nEllipseParams);
            else
                % Find the center of the pupil boundary points, place the boundary
                % points in a matrix and shift them to the center position
                xCenter=mean(Xc); yCenter=mean(Yc);
                centerMatrix = repmat([xCenter'; yCenter'], 1, length(Xc));
                
                % Rotate the data and split in half through the center
                pFitTransparentSplit=NaN(2,p.Results.nSplits,nEllipseParams);
                for ss=1:p.Results.nSplits
                    theta=((pi/2)/p.Results.nSplits)*ss;
                    forwardPoints = feval(returnRotMat,theta) * ([Xc,Yc]' - centerMatrix) + centerMatrix;
                    splitIdx1 = find((forwardPoints(1,:) < median(forwardPoints(1,:))))';
                    splitIdx2 = find((forwardPoints(1,:) >= median(forwardPoints(1,:))))';
                    pFitTransparentSplit(1,ss,:) = ...
                        feval(obtainPupilLikelihood,Xc(splitIdx1), Yc(splitIdx1));
                    pFitTransparentSplit(2,ss,:) = ...
                        feval(obtainPupilLikelihood,Xc(splitIdx2), Yc(splitIdx2));
                end % loop through splits
                
                % Calculate the SD of the parameters across splits, scaling by
                % sqrt(2) to roughly account for our use of just half the data
                pInitialFitSplitsSD=nanstd(reshape(pFitTransparentSplit,ss*2,nEllipseParams))/sqrt(2);
            end % check if we want to do splits
            
            % Obtain the SD of the parameters through a bootstrap resample
            % of data points if requested
            if p.Results.nBoots == 0
                pInitialFitBootsSD=NaN(1,nEllipseParams);
            else
                bootOptSet = statset('UseParallel',p.Results.useParallel);
                pInitialFitBootsSD = nanstd(bootstrp(p.Results.nBoots,obtainPupilLikelihood,Xc,Yc,'Options',bootOptSet));
            end % check if we want to do bootstraps
            
        end % check if there are pupil boundary data to be fit
        
        % store results
        loopVar_pInitialFitTransparent(ii,:) = pInitialFitTransparent';
        loopVar_pInitialFitHessianSD(ii,:) = pInitialFitHessianSD';
        loopVar_pInitialFitSplitsSD(ii,:) = pInitialFitSplitsSD';
        loopVar_pInitialFitBootsSD(ii,:) = pInitialFitBootsSD';
        loopVar_pInitialFitError(ii) = pInitialFitError;
    catch ME
        warning ('Error while processing frame: %d', ii)
        rethrow(ME)
    end % try catch
end % loop over frames

% gather the loop vars into the ellipse structure
pupilData.pInitialFitTransparent = loopVar_pInitialFitTransparent;
pupilData.pInitialFitHessianSD = loopVar_pInitialFitHessianSD;
pupilData.pInitialFitSplitsSD = loopVar_pInitialFitSplitsSD;
pupilData.pInitialFitBootsSD = loopVar_pInitialFitBootsSD;
pupilData.pInitialFitError = loopVar_pInitialFitError;

%% Clean up and save

% alert the user that we are done with the fit loop
if strcmp(p.Results.verbosity,'full')
    toc
    fprintf('\n');
end

% add a meta field with analysis details
pupilData.meta.fitPupilPerimeter = p.Results;
pupilData.meta.fitPupilPerimeter.coordinateSystem = 'intrinsicCoordinates(pixels)';

% save the ellipse fit results
save(p.Results.pupilFileName,'pupilData')


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


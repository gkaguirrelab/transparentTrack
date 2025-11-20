function [pupilData, priorPupilRadius] = smoothPupilTime(perimeterFileName, pupilFileName, sceneGeometryFileName, varargin)
% Constrain eye pose using a smoothly varying pupil radius
%
% Syntax:
%  [pupilData] = (newField)(perimeterFileName, pupilFileName, sceneGeometryFileName)
%
% Description:
%   Re-fits the eyePose values based upon a temporally smoothed version of
%   the prior estimate of pupil radius.
%
% Notes:
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
%                           perimeter data.
%   pupilFileName         - Full path to the .mat file that contains the
%                           pupil data to be smoothed. This file will be
%                           over-written by the output.
%   sceneGeometryFileName - Full path to the .mat file that contains the
%                           sceneGeometry to be used.
%
% Optional key/value pairs (display and I/O):
%  'verbose'              - Logical. Default false.
%
% Optional key/value pairs (flow control)
%  'nFrames'              - Analyze fewer than the total number of frames.
%  'useParallel'          - If set to true, use the Matlab parallel pool
%  'nWorkers'             - Specify the number of workers in the parallel
%                           pool. If undefined the default number will be
%                           used.
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
%  'glintFileName'        - Full path to a glint file. When available, the
%                           glint is used to constrain the eyePose that is
%                           found to fit the pupil perimeter.
%  'eyePoseLB'            - Lower bound on the eyePose
%  'eyePoseUB'            - Upper bound on the eyePose
%  'exponentialTauParam'  - The time constant (in video frames) of the
%                           decaying exponential weighting function for
%                           pupil radius.
%  'likelihoodErrorMultiplier' - The SD of the parameters estimated for
%                           each frame are computed as the product of the
%                           RMSE of the ellipse fits to the pupil perimeter
%                           points, the non-linear non-uniformity of the
%                           distribution of the points in space, and by
%                           this value. Typically set to ~4 to result in an
%                           SD of 1 when the fit of the points is good.
%                           Make this value larger to increase the
%                           influence of the prior, and smaller to increase
%                           the influence of the measurement.
%  'fitLabel'             - Identifies the field in pupilData that contains
%                           the ellipse fit params for which the search
%                           will be conducted.
%  'fixedPriorPupilRadius' - Scalar that provides the mean (in mm) of the
%                           expected radius of the pupil aperture during
%                           this acquisition. If set to empty (the default)
%                           the routine derives this values from the
%                           sceneConstrained fit results.
%  'relativeCameraPositionFileName' - Char. This is the full path to a
%                           relativeCameraPosition.mat file that provides
%                           the relative position of the camera at each
%                           video frame relative to the initial position of
%                           the camera.
%
% Outputs:
%   pupilData             - A structure with multiple fields corresponding
%                           to the parameters, SDs, and errors of the
%                           initial and final ellipse fits.
%

%% Parse vargin for options passed here
p = inputParser; p.KeepUnmatched = true; p.PartialMatching = false;

% Required
p.addRequired('perimeterFileName',@ischar);
p.addRequired('pupilFileName',@ischar);
p.addRequired('sceneGeometryFileName',@ischar);

% Optional display and I/O params
p.addParameter('verbose',false,@islogical);

% Optional flow control params
p.addParameter('startFrame',1,@isnumeric);
p.addParameter('useParallel',false,@islogical);
p.addParameter('nWorkers',[],@(x)(isempty(x) | isnumeric(x)));

% Optional environment parameters
p.addParameter('tbSnapshot',[],@(x)(isempty(x) | isstruct(x)));
p.addParameter('timestamp',char(datetime('now')),@ischar);
p.addParameter('hostname',char(java.lang.System.getProperty('user.name')),@ischar);
p.addParameter('username',char(java.net.InetAddress.getLocalHost.getHostName),@ischar);

% Optional fitting params
p.addParameter('glintFileName',[],@(x)(isempty(x) || ischar(x)));
p.addParameter('eyePoseLB',[-89,-89,0,0.1],@isnumeric);
p.addParameter('eyePoseUB',[89,89,0,5],@isnumeric);
p.addParameter('cameraTransBounds',[0; 0; 0],@isnumeric);
p.addParameter('glintTol',1.5,@isnumeric);
p.addParameter('priorPupilRadius',[],@isnumeric);
p.addParameter('priorPupilBound',0.125,@isnumeric);
p.addParameter('confidenceThreshold',0.75,@isnumeric);
p.addParameter('pupilVarThresh',0.05,@isnumeric);
p.addParameter('ellipseRMSEThresh',1,@isnumeric);
p.addParameter('minPerimPoints',5,@isnumeric);
p.addParameter('justReturnPrior',false,@islogical);
p.addParameter('currField','',@ischar);


%% Parse and check the parameters
p.parse(perimeterFileName, pupilFileName, sceneGeometryFileName, varargin{:});

nEllipseParams=5; % 5 params in the transparent ellipse form
nEyePoseParams=4; % 4 eyePose values (azimuth, elevation, torsion, radius)
radiusIdx = 4; % The 4th eyePose entry holds the radius value
nHeadTransParams=3; % [horizontal; vertical; depth]
cameraTransBounds = [0;0;0];
fps = 120;

%% Load and check data
% Load the pupil perimeter data
load(perimeterFileName,'perimeter');

% Load the pupil data
load(pupilFileName,'pupilData');

% Load the sceneGeometry file
% Need to do the inelegant loading approach to keep parloop happy
dataLoad = load(sceneGeometryFileName);
sceneGeometry = dataLoad.sceneGeometry;
clear dataLoad

% Load the glint file if passed
glintTol = p.Results.glintTol;
if ~isempty(p.Results.glintFileName)
    load(p.Results.glintFileName,'glintData');
    % Sometimes the glint file has fewer frames than the perimeter data. If
    % so, make sure that we discard any perimeter frames that do not have a
    % corresponding glint
    nGlintFrames = size(glintData.X,1);
    if nGlintFrames < size(perimeter.data,1)
        perimeter.data = perimeter.data(1:nGlintFrames);
    end
else
    glintData = [];
end

% determine how many frames we will process
startFrame = p.Results.startFrame;
nFrames=size(perimeter.data,1);

% determine the current field of the pupilData structure
if ~isempty(p.Results.currField)
    currField = p.Results.currField;
else
    currField = pupilData.currentField;
end

%% Derive a smooth version of the pupil size from a previous analysis
if isempty(p.Results.priorPupilRadius)
    y = pupilData.(currField).eyePoses.values(:,4);
    rmse = pupilData.(currField).ellipses.RMSE;

    % Measure the moving standard deviation of pupil size
    pupilMovSD = movstd(y,fps/4,"omitnan");

    % Identify the "good" pupil size values
    goodPupilTimeIdx = and(pupilMovSD < p.Results.pupilVarThresh,rmse < p.Results.ellipseRMSEThresh);

    % A version of the pupil radius, with nans for the not-good values
    yMod = y; yMod(~goodPupilTimeIdx)=nan;

    % Replace nan values at the start and end with the median pupil size
    % across the vector
    firstNonNan = find(~isnan(yMod),1);
    if firstNonNan>1
        yMod(1:firstNonNan-1) = median(yMod,"omitmissing");
    end
    lastNonNan = find(~isnan(yMod),1,'last');
    if lastNonNan < length(yMod)
        yMod(lastNonNan+1:end) = median(yMod,"omitmissing");
    end

    % Fill missing values with linear interpolation
    yMod = fillmissing(yMod,"linear");

    % Smooth the final vector
    priorPupilRadius = smoothdata(yMod,"loess",fps);
else
    priorPupilRadius = p.Results.priorPupilRadius;
end

% Check if we just wanted to examine the prior for these settings
if p.Results.justReturnPrior
    return
end

% Derive a bound based upon the SD of the prior (if not otherwise defined)
if isempty(p.Results.priorPupilBound)
    priorPupilBound = std(priorPupilRadius);
else
    priorPupilBound = p.Results.priorPupilBound;
end


%% Set up the parallel pool
if p.Results.useParallel
    nWorkers = startParpool( p.Results.nWorkers, p.Results.verbose );
else
    nWorkers=0;
end

% Recast perimeter.data into a sliced cell array to reduce parfor
% broadcast overhead
frameCellArray = perimeter.data;
clear perimeter

% Set-up other variables to be non-broadcast
verbose = p.Results.verbose;
eyePoseLB = p.Results.eyePoseLB;
eyePoseUB = p.Results.eyePoseUB;
confidenceThreshold = p.Results.confidenceThreshold;
minPerimPoints = p.Results.minPerimPoints;

%% Perform the calculation across frames

% Alert the user
if p.Results.verbose
    tic
    fprintf(['Bayesian smoothing. Started ' char(datetime('now')) '\n']);
    fprintf('| 0                      50                   100%% |\n');
    fprintf('.\n');
end

% Store the warning state
warnState = warning();

% Loop through the frames
parfor (ii = startFrame:nFrames, nWorkers)
    %for ii = startFrame:nFrames

    % update progress
    if verbose
        if mod(ii,round((nFrames-startFrame)/50))==0
            fprintf('\b.\n');
        end
    end

    % initialize some variables so that their use is transparent to the
    % parfor loop
    ellipseParams = NaN(1,nEllipseParams);
    eyePose = NaN(1,nEyePoseParams);
    cameraTrans = zeros(nHeadTransParams,1);
    fitAtBound = false;
    fVal = nan;

    % get the boundary points
    Xp = frameCellArray{ii}.Xp;
    Yp = frameCellArray{ii}.Yp;

    % Check if we have a confidence field. If so, filter the points
    if isfield(frameCellArray{ii},'confidence')
        conf = frameCellArray{ii}.confidence;
        goodIdx = conf > confidenceThreshold;
        Xp = Xp(goodIdx);
        Yp = Yp(goodIdx);
    end

    % if this frame has data, proceed to fit
    if ~isempty(Xp) &&  ~isempty(Yp)

        % Require a minimum number of valid perimeter points
        if length(Xp) > minPerimPoints

            % Re-fit the ellipse with the radius constrained to the
            % posterior value. Pass the prior azimuth and elevation as x0.
            lb_pin = eyePoseLB;
            ub_pin = eyePoseUB;
            lb_pin(radiusIdx) = priorPupilRadius(ii)-priorPupilBound;
            ub_pin(radiusIdx) = priorPupilRadius(ii)+priorPupilBound;
            x0 = pupilData.(currField).eyePoses.values(ii,:);
            x0(radiusIdx)=priorPupilRadius(ii);

            % If the prior fit was not good, set the x0 guess for azi and
            % ele to the mean of the bounds. This avoids getting stuck in
            % local bad minima.
            if ~goodPupilTimeIdx(ii)
                x0(1) = mean([eyePoseLB(1); eyePoseUB(1)]);
                x0(2) = mean([eyePoseLB(2); eyePoseUB(2)]);
            end

            % If we have any nans in the x0 guess, clear it
            if any(isnan(x0))
                x0 = [];
            end

            % If we have glintData, extract the glintCoord, and allow
            % non-zero bounds on the cameraTrans search. If no glint data,
            % then lock the cameraTransBounds to zero.
            if ~isempty(glintData)
                glintCoord = [glintData.X(ii,:), glintData.Y(ii,:)];
            else
                glintCoord = [];
            end

            if isempty(glintCoord) || any(isnan(glintCoord))
                thisFrameCameraTransBounds = [0; 0; 0];
            else
                thisFrameCameraTransBounds = cameraTransBounds;
            end

            % Turn off warnings that can arise when fitting bad frames
            warning('off','projectModelEye:rayTracingError');
            warning('off','projectModelEye:ellipseFitFailed');
            warning('off','gkaModelEye:pupilEllipseFit');

            % Perform the fit
            [eyePose, cameraTrans, fVal, ellipseParams, fitAtBound] = ...
                eyePoseEllipseFit(Xp, Yp, glintCoord, sceneGeometry, ...
                'glintTol',glintTol,'cameraTransX0',cameraTrans, ...
                'cameraTransBounds',thisFrameCameraTransBounds, ...
                'eyePoseLB', lb_pin, 'eyePoseUB', ub_pin, 'eyePoseX0', x0);

            % Restore the warning state
            warning(warnState);

        end % Check that we have enough perimeter points

    end % check if there are any perimeter points to fit

    % store results
    loopVar_eyePose(ii,:) = eyePose;
    loopVar_cameraTrans(ii,:) = cameraTrans;
    loopVar_fVal(ii) = fVal;
    loopVar_ellipseParams(ii,:) = ellipseParams;
    loopVar_fitAtBound(ii) = fitAtBound;

end % loop over frames to calculate the posterior

% report completion of Bayesian analysis
if p.Results.verbose
    toc
    fprintf('\n');
end

%% Clean up and save the fit results
if contains(currField,'smoothPupilTime')
    parts = split(currField,'_');
    newField = sprintf('smoothPupilTime_%02d',str2double(parts{2})+1);
else
    newField = 'smoothPupilTime_01';
end

% Clear out any prior results in the (newField) field
pupilData.(newField) = [];

% gather the loop vars into the ellipses field
pupilData.(newField).ellipses.values=loopVar_ellipseParams;
pupilData.(newField).ellipses.RMSE=loopVar_fVal';
pupilData.(newField).ellipses.meta.ellipseForm = 'transparent';
pupilData.(newField).ellipses.meta.labels = {'x','y','area','eccentricity','theta'};
pupilData.(newField).ellipses.meta.units = {'pixels','pixels','squared pixels','non-linear eccentricity','rads'};
pupilData.(newField).ellipses.meta.coordinateSystem = 'intrinsic image';

% gather the loop vars into the eyePoses field
pupilData.(newField).eyePoses.priorPupilRadius=priorPupilRadius;
pupilData.(newField).eyePoses.values=loopVar_eyePose;
pupilData.(newField).eyePoses.fitAtBound = loopVar_fitAtBound';
pupilData.(newField).eyePoses.meta.labels = {'azimuth','elevation','torsion','pupil radius'};
pupilData.(newField).eyePoses.meta.units = {'deg','deg','deg','mm'};
pupilData.(newField).eyePoses.meta.coordinateSystem = 'head fixed (extrinsic)';

% add a meta field with analysis details
pupilData.(newField).meta = p.Results;

% Update the relativeCameraPosition
relativeCameraPosition.(newField).values = loopVar_cameraTrans';
relativeCameraPosition.(newField).meta = p.Results;
relativeCameraPosition.(newField) = '(newField)';

% Store the identity of the most recently produced field of data
pupilData.currentField = newField;

% save the pupilData
save(pupilFileName,'pupilData')

end % function

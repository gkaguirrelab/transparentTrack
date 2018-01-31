function sceneGeometry = estimateSceneGeometry(pupilFileName, sceneGeometryFileName, varargin)
% Estimate eye radius and eye center given a set of image plane ellipses
%
% Description:
%   This function searches over a set of ellipses from the passed pupil
%   file(s) to estimate the extrinsic camera translation vector. The search
%   attempts to minimize the error associated with the prediction of the
%   center of ellipses in the image plane while constraining the shape of
%   these ellipses to match that predicted by the projection.
%
%   We are aware that the center of an ellipse on the image plane is not at
%   the same location as the projection of the center of the pupil on the
%   image plane (see, for example,  Ahn, Sung Joon, H. J. Warnecke, and
%   Rüdiger Kotowski. "Systematic geometric image measurement errors of
%   circular object targets: Mathematical formulation and correction." The
%   Photogrammetric Record 16.93 (1999): 485-502.). The modeling solution
%   implemented here accounts for this property, as we implement a full,
%   numeric forward projection of the pupil circle to the image plane,
%   including accounting for refraction at the cornea.
%
% Components of the projection model:
%   intrinsicCameraMatrix - This matrix has the form:
%
%       [fx  s x0]
%       [0  fy y0]
%       [0   0  1]
%
%   where fx and fy are the focal lengths of the camera in the x and y
%   image dimensions, s is the axis skew, and x0, y0 define the principle
%   offset point. For a camera sensor with square pixels, fx = fy. Ideally,
%   skew should be zero. The principle offset point should be in the center
%   of the sensor. Units are traditionally in pixels. These values can be
%   empirically measured for a camera using a calibration approach
%   (https://www.mathworks.com/help/vision/ref/cameramatrix.html).
%   The default values are those that we measured by calibration of the
%   LiveTrack AV, 12M-i camera.
%
%   radialDistortionVector - A two element vector of the form:
%
%       [k1 k2]
%
%   that models the radial distortion introduced the lens. This is an
%   empirically measured property of the camera system and so these
%   parameters are locked here. The default values are those that we
%   measured by calibration of the LiveTrack AV, 12M-i camera.
%
%   extrinsicTranslationVector - a vector of the form:
%
%       [x]
%       [y]
%       [z]
%
%   with the values specifying the location (horizontal, vertical, and
%   depth, respectively) of the principle offset point of the camera in mm
%   relative to the scene coordinate system. We define the origin of the
%   scene coordinate system to be x=0, y=0 along the optical axis of the
%   eye, and z=0 to be the apex of the corneal surface. The default values
%   and bounds can be guided by knowledge of the physical arrangement of
%   the subject and recording apparatus. These values are the target of the
%   search.
%
%   extrinsicRotationMatrix - A 3x3 matrix with the fixed values of:
%
%       [1  0  0]
%       [0 -1  0]
%       [0  0 -1]
%
%   These values invert the axes of the scene coordinate system to produce
%   the direction conventions of the image coordinate system. As the
%   projection of pupil circles from scene to image is invariant to
%   rotations of the camera matrix, these values are locked.
%
%   primaryPosition - A 1x3 vector of:
%
%       [eyeAzimuth, eyeElevation, eyeTorsion]
%
%   that specifies the rotation angles (in head fixed axes) for which the
%   eye is in primary position (as defined by Listing's Law). These are set
%   to zero and are not modified in this routine.
%
%   constraintTolerance - A scalar. The inverse projection from ellipse on
%   the image plane to eye params (azimuth, elevation) imposes a constraint
%   on how well the solution must match the shape of the ellipse (defined
%   by ellipse eccentricity and theta) and the area of the ellipse. This
%   constraint is expressed as a proportion of error, relative to either
%   the largest possible area in ellipse shape or an error in area equal to
%   the area of the ellipse itself (i.e., unity). If the constraint is made
%   too stringent, then in an effort to perfectly match the shape of the
%   ellipse, error will increase in matching the position of the ellipse
%   center. It should be noted that even in a noise-free simulation, it is
%   not possible to perfectly match ellipse shape and area while matching
%   ellipse center position, as the shape of the projection of the pupil
%   upon the image plane deviates from perfectly elliptical due to
%   perspective effects. We find that a value in the range 0.01 - 0.03
%   provides an acceptable compromise in empirical data.
%
%   eye -  This is itself a structure that is returned by the function
%   modelEyeParameters(). The parameters define the anatomical properties
%   of the eye, including the size and shape of the anterior and posterior
%   chamber. These parameters are adjusted for the measured spherical
%   refractive error of the subject and (optionally) measured axial length.
%
% Inputs:
%	pupilFileName         - Full path to a pupilData file, a cell array
%                           of such paths, or a pupilData structure itself.
%                           If a single path, the pupilData file is loaded.
%                           If a cell array, the ellipse data from each
%                           pupilData file is loaded and concatenated. If
%                           set to empty, a sceneGeometry structure with
%                           default values will be returned.
%   sceneGeometryFileName - Full path to the file in which the
%                           sceneGeometry data should be saved
%
% Optional key/value pairs (display and I/O):
%  'verbosity'            - Level of verbosity. [none, full]
%  'sceneDiagnosticPlotFileName' - Full path (including suffix) to the
%                           location where a diagnostic plot of the
%                           sceneGeometry calculation is to be saved. If
%                           left empty, then no plot will be saved.
%
% Optional key/value pairs (flow control)
%  'useParallel'          - If set to true, use the Matlab parallel pool
%  'nWorkers'             - Specify the number of workers in the parallel
%                           pool. If undefined the default number will be
%                           used.
%  'tbtbProjectName'      - The workers in the parallel pool are configured
%                           by issuing a tbUseProject command for the
%                           project specified here.
%
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
% Optional key/value pairs (analysis)
%  'radialDistortionVector' - 1x2 vector of radial distortion parameters
%  'intrinsicCameraMatrix' - 3x3 matrix
%  'extrinsicRotationMatrix' - 3x3 matrix
%  'extrinsicTranslationVector' - 3x1 vector
%  'extrinsicTranslationVectorLB' - 3x1 vector
%  'extrinsicTranslationVectorUB' - 3x1 vector
%  'spectacleRefractionDiopters' - Scalar. A negative value is appropriate
%                           for a myopic subject.
%  'axialLength'          - Scalar, defaults to empty. When set, this
%                           causes the posterior chamber of the model eye
%                           (and rotation center) to be adjusted to match
%                           the passed empirical axial length (in mm).
%  'primaryPosition'      - 1x3 vector
%  'constraintTolerance'  - Scalar. Range 0-1. Typical value 0.01 - 0.03
%  'eyeParamsLB/UB'       - Upper and lower bounds on the eyeParams
%                           [azimuth, elevation, torsion, pupil radius].
%                           The torsion value is unusued and is bounded to
%                           zero. Biological limits in eye rotation and
%                           pupil size would suggest boundaries of [±35,
%                           ±25, 0, 0.25-5]. Note, however, that these
%                           angles are relative to the center of
%                           projection, not the primary position of the
%                           eye. Therefore, in circumstances in which the
%                           camera is viewing the eye from an off-center
%                           angle, the bounds will need to be shifted
%                           accordingly.
%  'fitLabel'             - Identifies the field in pupilData that contains
%                           the ellipse fit params for which the search
%                           will be conducted.
%  'ellipseArrayList'     - A vector of frame numbers (indexed from 1)
%                           which identify the llipses to be used for the
%                           estimation of scene geometry. If left empty,
%                           a list of ellipses will be generated.
%  'nBinsPerDimension'    - Scalar. Defines the number of divisions with
%                           which the ellipse centers are binned.
%
% Outputs
%	sceneGeometry         - A structure that contains the components of the
%                           projection model.
%


%% input parser
p = inputParser; p.KeepUnmatched = true;

% Required
p.addRequired('pupilFileName',@(x)(isempty(x) | isstruct(x) | iscell(x) | ischar(x)));
p.addRequired('sceneGeometryFileName',@(x)(isempty(x) | ischar(x)));

% Optional display and I/O params
p.addParameter('verbosity', 'none', @isstr);
p.addParameter('sceneDiagnosticPlotFileName', '', @(x)(isempty(x) | ischar(x)));

% Optional flow control params
p.addParameter('useParallel',false,@islogical);
p.addParameter('nWorkers',[],@(x)(isempty(x) | isnumeric(x)));
p.addParameter('tbtbRepoName','transparentTrack',@ischar);

% Optional environment parameters
p.addParameter('tbSnapshot',[],@(x)(isempty(x) | isstruct(x)));
p.addParameter('timestamp',char(datetime('now')),@ischar);
p.addParameter('username',char(java.lang.System.getProperty('user.name')),@ischar);
p.addParameter('hostname',char(java.net.InetAddress.getLocalHost.getHostName),@ischar);

% Optional analysis params
p.addParameter('radialDistortionVector',[-0.3517 3.5353],@isnumeric);
p.addParameter('intrinsicCameraMatrix',[2627.0 0 338.1; 0 2628.1 246.2; 0 0 1],@isnumeric);
p.addParameter('extrinsicRotationMatrix',[1 0 0; 0 -1 0; 0 0 -1],@isnumeric);
p.addParameter('extrinsicTranslationVector',[0; 0; 120],@isnumeric);
p.addParameter('extrinsicTranslationVectorLB',[-10; -10; 100],@isnumeric);
p.addParameter('extrinsicTranslationVectorUB',[10; 10; 180],@isnumeric);
p.addParameter('spectacleRefractionDiopters',0,@isnumeric);
p.addParameter('axialLength',[],@(x)(isempty(x) | isnumeric(x)));
p.addParameter('eyeLaterality','right',@ischar);
p.addParameter('primaryPosition',[0 0 0],@isnumeric);
p.addParameter('constraintTolerance',0.02,@isnumeric);
p.addParameter('eyeParamsLB',[-35,-25,0,0.25],@(x)(isempty(x) | isnumeric(x)));
p.addParameter('eyeParamsUB',[35,25,0,4],@(x)(isempty(x) | isnumeric(x)));
p.addParameter('fitLabel','initial',@ischar);
p.addParameter('ellipseArrayList',[],@(x)(isempty(x) | isnumeric(x)));
p.addParameter('nBinsPerDimension',10,@isnumeric);

% parse
p.parse(pupilFileName, sceneGeometryFileName, varargin{:})


%% Announce we are starting
if strcmp(p.Results.verbosity,'full')
    tic
    fprintf(['Estimating scene geometry from pupil ellipses. Started ' char(datetime('now')) '\n']);
end

%% Create initial sceneGeometry structure and ray tracing functions
% Announce
if strcmp(p.Results.verbosity,'full')
    fprintf('Creating initial scene geometry.\n');
end
initialSceneGeometry.radialDistortionVector = p.Results.radialDistortionVector;
initialSceneGeometry.intrinsicCameraMatrix = p.Results.intrinsicCameraMatrix;
initialSceneGeometry.extrinsicTranslationVector = p.Results.extrinsicTranslationVector;
initialSceneGeometry.extrinsicRotationMatrix = p.Results.extrinsicRotationMatrix;
initialSceneGeometry.primaryPosition = p.Results.primaryPosition;
initialSceneGeometry.constraintTolerance = p.Results.constraintTolerance;
initialSceneGeometry.eye = modelEyeParameters(...
    'spectacleRefractionDiopters',p.Results.spectacleRefractionDiopters,...
    'axialLength',p.Results.axialLength,...
    'eyeLaterality',p.Results.eyeLaterality);

% Return the initialSceneGeometry if pupilFileName is empty
if isempty(pupilFileName)
    sceneGeometry = initialSceneGeometry;
    return
end

% Assemble the ray tracing functions
if strcmp(p.Results.verbosity,'full')
    fprintf('Assembling ray tracing functions.\n');
end
[rayTraceFuncs] = assembleRayTraceFuncs( initialSceneGeometry );


%% Set up the parallel pool
if p.Results.useParallel
    % If a parallel pool does not exist, attempt to create one
    poolObj = gcp('nocreate');
    if isempty(poolObj)
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
        nWorkers = poolObj.NumWorkers;
    end
else
    nWorkers=0;
end


%% Load pupil data
if iscell(pupilFileName)
    ellipses = [];
    ellipseFitSEM = [];
    for cc = 1:length(pupilFileName)
        load(pupilFileName{cc})
        ellipses = [ellipses;pupilData.(p.Results.fitLabel).ellipses.values];
        ellipseFitSEM = [ellipseFitSEM; pupilData.(p.Results.fitLabel).ellipses.RMSE];
    end
end
if ischar(pupilFileName)
    load(pupilFileName)
    ellipses = pupilData.(p.Results.fitLabel).ellipses.values;
    ellipseFitSEM = pupilData.(p.Results.fitLabel).ellipses.RMSE;
end
if isstruct(pupilFileName)
    pupilData = pupilFileName;
    ellipses = pupilData.(p.Results.fitLabel).ellipses.values;
    ellipseFitSEM = pupilData.(p.Results.fitLabel).ellipses.RMSE;
end


%% Identify the ellipses that will guide the sceneGeometry estimation
% If not supplied, we will generate a list of ellipses to use for the
% estimation.
if ~isempty(p.Results.ellipseArrayList)
    ellipseArrayList = p.Results.ellipseArrayList;
    Xedges = [];
    Yedges = [];
else
    if strcmp(p.Results.verbosity,'full')
        fprintf('Selecting ellipses to guide the search.\n');
    end
    
    % First we divide the ellipse centers amongst a set of 2D bins across image
    % space. We will ultimately minimize the fitting error across bins
    [ellipseCenterCounts,Xedges,Yedges,binXidx,binYidx] = ...
        histcounts2(ellipses(:,1),ellipses(:,2),p.Results.nBinsPerDimension);
    
    % Anonymous functions for row and column identity given array position
    rowIdx = @(b) fix( (b-1) ./ (size(ellipseCenterCounts,2)) ) +1;
    colIdx = @(b) 1+mod(b-1,size(ellipseCenterCounts,2));
    
    % Create a cell array of index positions corresponding to each of the 2D
    % bins
    idxByBinPosition = ...
        arrayfun(@(b) find( (binXidx==rowIdx(b)) .* (binYidx==colIdx(b)) ),1:1:numel(ellipseCenterCounts),'UniformOutput',false);
    
    % Identify which bins are not empty
    filledBinIdx = find(~cellfun(@isempty, idxByBinPosition));
    
    % Identify the ellipses in each filled bin with the lowest fit SEM
    [~, idxMinErrorEllipseWithinBin] = arrayfun(@(x) nanmin(ellipseFitSEM(idxByBinPosition{x})), filledBinIdx, 'UniformOutput', false);
    returnTheMin = @(binContents, x)  binContents(idxMinErrorEllipseWithinBin{x});
    ellipseArrayList = cellfun(@(x) returnTheMin(idxByBinPosition{filledBinIdx(x)},x),num2cell(1:1:length(filledBinIdx)));
end


%% Generate the errorWeights
errorWeights= ellipseFitSEM(ellipseArrayList);
errorWeights = 1./errorWeights;
errorWeights=errorWeights./mean(errorWeights);


%% Perform the search
if strcmp(p.Results.verbosity,'full')
    fprintf('Searching over camera translation and eye rotation center.\n');
end
% Call out to the local function that performs the serach
sceneGeometry = ...
    performSceneSearch(initialSceneGeometry, rayTraceFuncs, ...
    ellipses(ellipseArrayList,:), ...
    errorWeights, ...
    p.Results.extrinsicTranslationVectorLB, ...
    p.Results.extrinsicTranslationVectorUB, ...
    p.Results.eyeParamsLB, ...
    p.Results.eyeParamsUB, ...
    nWorkers);

% add additional search and meta field info to sceneGeometry
sceneGeometry.meta.estimateGeometry.parameters = p.Results;
sceneGeometry.meta.estimateGeometry.search.ellipseArrayList = ellipseArrayList';


%% Save the sceneGeometry file
if ~isempty(sceneGeometryFileName)
    save(sceneGeometryFileName,'sceneGeometry');
end


%% Create a sceneGeometry plot
if ~isempty(p.Results.sceneDiagnosticPlotFileName)
    if strcmp(p.Results.verbosity,'full')
        fprintf('Creating a sceneGeometry diagnostic plot.\n');
    end
    saveSceneDiagnosticPlot(...
        ellipses(ellipseArrayList,:),...
        Xedges, Yedges,...
        p.Results.eyeParamsLB, ...
        p.Results.eyeParamsUB, ...
        sceneGeometry,...
        rayTraceFuncs,...
        p.Results.sceneDiagnosticPlotFileName)
end


%% alert the user that we are done with the routine
if strcmp(p.Results.verbosity,'full')
    toc
    fprintf('\n');
end


end % main function



%% LOCAL FUNCTIONS

function sceneGeometry = performSceneSearch(initialSceneGeometry, rayTraceFuncs, ellipses, errorWeights, translationVectorLB, translationVectorUB, eyeParamsLB, eyeParamsUB, nWorkers)
% Pattern search for best fitting sceneGeometry parameters
%
% Description:
%   The routine searches for parameters of the extrinsic translation vector
%   of the camera that best models the locations of the centers of ellipses
%   found on the image plane, given the constraint that the ellipse shape
%   (and area) must also match the prediction of the forward model. The
%   passed sceneGeometry structure is used as the starting point for the
%   search. Across each iteration of the search, a candidate sceneGeometry
%   is assembled from the current values of the parameters. This
%   sceneGeometry is then used in the inverse pupil projection model. The
%   inverse projection searches for an eye azimuth, elevation, and pupil
%   radius that, given the sceneGeometry, best accounts for the parameters
%   of the target ellipse on the image plane. This inverse search attempts
%   to minimize the distance bewteen the centers of the predicted and
%   targeted ellipse on the image plane, while satisfying non-linear
%   constraints upon matching the shape (eccentricity and theta) and area
%   of the ellipses. Only when the translation vector is correctly
%   specified will the inverse pupil projection model be able to
%   simultaneouslty match the center and shape of the ellipse on the image
%   plane.
%
%   The iterative search across sceneGeometry parameters attempts to
%   minimize the L(norm) of the distances between the targeted and modeled
%   centers of the ellipses. In the calculation of this objective functon,
%   each distance error is weighted. The error weight is derived from the
%   accuracy with which the boundary points of the pupil in the image plane
%   are fit by an unconstrained ellipse.
%
%   We find that the patternsearch optimization gives the best results in
%   the shortest period.
%

% Set the error form
errorForm = 'RMSE';

% Extract the initial search point from initialSceneGeometry
x0 = initialSceneGeometry.extrinsicTranslationVector;

% Define search options
options = optimoptions(@patternsearch, ...
    'Display','iter',...
    'AccelerateMesh',false,...
    'Cache','on',...
    'CompleteSearch','on',...
    'FunctionTolerance',1e-6);

% Define anonymous functions for the objective and constraint
objectiveFun = @objfun; % the objective function, nested below

% Define nested variables for within the search
centerDistanceErrorByEllipse=zeros(size(ellipses,1),1);
shapeErrorByEllipse=zeros(size(ellipses,1),1);
areaErrorByEllipse=zeros(size(ellipses,1),1);

[x, fVal] = patternsearch(objectiveFun, x0,[],[],[],[],fliplr(translationVectorLB),fliplr(translationVectorUB),[],options);
% Nested function computes the objective for the patternsearch
    function fval = objfun(x)
        % Assemble a candidate sceneGeometry structure
        candidateSceneGeometry = initialSceneGeometry;
        candidateSceneGeometry.extrinsicTranslationVector = x;
        % For each ellipse, perform the inverse projection from the ellipse
        % on the image plane to eyeParams. We retain the errors from the
        % inverse projection and use these to assemble the objective
        % function. We parallelize the computation across ellipses.
        parfor (ii = 1:size(ellipses,1), nWorkers)
            [~, ~, centerDistanceErrorByEllipse(ii), shapeErrorByEllipse(ii), areaErrorByEllipse(ii)] = ...
                pupilProjection_inv(...
                ellipses(ii,:),...
                candidateSceneGeometry, rayTraceFuncs, ...
                'eyeParamsLB',eyeParamsLB,...
                'eyeParamsUB',eyeParamsUB...
                );
        end
        % Now compute objective function as the RMSE of the distance
        % between the taget and modeled ellipses
        switch errorForm
            case 'SSE'
                fval=sum((centerDistanceErrorByEllipse.*(shapeErrorByEllipse.*100+1).*(areaErrorByEllipse.*100+1).*errorWeights).^2);
            case 'RMSE'
                fval = mean((centerDistanceErrorByEllipse.*(shapeErrorByEllipse.*100+1).*(areaErrorByEllipse.*100+1).*errorWeights).^2).^(1/2);
            otherwise
                error('I do not recognize that error form');
        end
        
    end


% Assemble the sceneGeometry file to return
sceneGeometry.radialDistortionVector = initialSceneGeometry.radialDistortionVector;
sceneGeometry.intrinsicCameraMatrix = initialSceneGeometry.intrinsicCameraMatrix;
sceneGeometry.extrinsicTranslationVector = x;
sceneGeometry.extrinsicRotationMatrix = initialSceneGeometry.extrinsicRotationMatrix;
sceneGeometry.primaryPosition = initialSceneGeometry.primaryPosition;
sceneGeometry.constraintTolerance = initialSceneGeometry.constraintTolerance;
sceneGeometry.eye = initialSceneGeometry.eye;
sceneGeometry.meta.estimateGeometry.search.options = options;
sceneGeometry.meta.estimateGeometry.search.errorForm = errorForm;
sceneGeometry.meta.estimateGeometry.search.initialSceneGeometry = initialSceneGeometry;
sceneGeometry.meta.estimateGeometry.search.ellipses = ellipses;
sceneGeometry.meta.estimateGeometry.search.errorWeights = errorWeights;
sceneGeometry.meta.estimateGeometry.search.sceneParamsLB = translationVectorLB;
sceneGeometry.meta.estimateGeometry.search.sceneParamsUB = translationVectorUB;
sceneGeometry.meta.estimateGeometry.search.eyeParamsLB = eyeParamsLB;
sceneGeometry.meta.estimateGeometry.search.eyeParamsUB = eyeParamsUB;
sceneGeometry.meta.estimateGeometry.search.fVal = fVal;
sceneGeometry.meta.estimateGeometry.search.centerDistanceErrorByEllipse = centerDistanceErrorByEllipse;
sceneGeometry.meta.estimateGeometry.search.shapeErrorByEllipse = shapeErrorByEllipse;
sceneGeometry.meta.estimateGeometry.search.areaErrorByEllipse = areaErrorByEllipse;

end % local search function


function [] = saveSceneDiagnosticPlot(ellipses, Xedges, Yedges, eyeParamsLB, eyeParamsUB, sceneGeometry, rayTraceFuncs, sceneDiagnosticPlotFileName)
% Creates and saves a plot that illustrates the sceneGeometry results
%
% Inputs:
%   ellipses              - An n x p array containing the p parameters of
%                           the n ellipses used to derive sceneGeometry
%   Xedges                - The X-dimension edges of the bins used to
%                           divide and select ellipses across the image.
%   Yedges                - The Y-dimension edges of the bins used to
%                           divide and select ellipses across the image.
%   eyeParamsLB, eyeParamsUB - Bounds for the eye params to be passed to
%                           pupilProjection_inv.
%   sceneGeometry         - The sceneGeometry structure
%   sceneDiagnosticPlotFileName - The full path (including .pdf suffix)
%                           to the location to save the diagnostic plot
%
% Outputs:
%   none
%

figHandle=figure('visible','off');
set(gcf,'PaperOrientation','landscape');

set(figHandle, 'Units','inches')
height = 6;
width = 11;

% the last two parameters of 'Position' define the figure size
set(figHandle, 'Position',[25 5 width height],...
    'PaperSize',[width height],...
    'PaperPositionMode','auto',...
    'Color','w',...
    'Renderer','painters'...     %recommended if there are no alphamaps
    );

%% Left panel -- distance error
subplot(3,3,[1 4]);

if ~isempty(Xedges)
    % plot the 2D histogram grid
    for xx = 1: length(Xedges)
        if xx==1
            hold on
        end
        plot([Xedges(xx) Xedges(xx)], [Yedges(1) Yedges(end)], '-', 'Color', [0.9 0.9 0.9], 'LineWidth', 0.5 );
    end
    for yy=1: length(Yedges)
        plot([Xedges(1) Xedges(end)], [Yedges(yy) Yedges(yy)], '-', 'Color', [0.9 0.9 0.9], 'LineWidth', 0.5);
    end
    binSpaceX = Xedges(2)-Xedges(1);
    binSpaceY = Yedges(2)-Yedges(1);
end

% plot the ellipse centers
scatter(ellipses(:,1),ellipses(:,2),'o','filled', ...
    'MarkerFaceAlpha',2/8,'MarkerFaceColor',[0 0 0]);
hold on

% get the predicted ellipse centers
[~, projectedEllipses] = ...
    arrayfun(@(x) pupilProjection_inv...
    (...
    ellipses(x,:),...
    sceneGeometry,...
    rayTraceFuncs,...
    'eyeParamsLB',eyeParamsLB,'eyeParamsUB',eyeParamsUB),...
    1:1:size(ellipses,1),'UniformOutput',false);
projectedEllipses=vertcat(projectedEllipses{:});

% plot the projected ellipse centers
scatter(projectedEllipses(:,1),projectedEllipses(:,2),'o','filled', ...
    'MarkerFaceAlpha',2/8,'MarkerFaceColor',[0 0 1]);

% connect the centers with lines
errorWeightVec=sceneGeometry.meta.estimateGeometry.search.errorWeights;
for ii=1:size(ellipses,1)
    lineAlpha = errorWeightVec(ii)/max(errorWeightVec);
    lineWeight = 0.5 + (errorWeightVec(ii)/max(errorWeightVec));
    ph=plot([projectedEllipses(ii,1) ellipses(ii,1)], ...
        [projectedEllipses(ii,2) ellipses(ii,2)], ...
        '-','Color',[1 0 0],'LineWidth', lineWeight);
    ph.Color(4) = lineAlpha;
end

% plot the estimated center of rotation of the eye
rotationCenterEllipse = pupilProjection_fwd([0 0 0 2], sceneGeometry, rayTraceFuncs);
plot(rotationCenterEllipse(1),rotationCenterEllipse(2), '+g', 'MarkerSize', 5);

% Calculate the plot limits
if ~isempty(Xedges)
    xPlotBounds = [Xedges(1)-binSpaceX Xedges(end)+binSpaceX];
    yPlotBounds = [Yedges(1)-binSpaceY Yedges(end)+binSpaceY];
else
    minX = min([projectedEllipses(:,1);ellipses(:,1)]);
    maxX = max([projectedEllipses(:,1);ellipses(:,1)]);
    minY = min([projectedEllipses(:,2);ellipses(:,2)]);
    maxY = max([projectedEllipses(:,2);ellipses(:,2)]);
    xPlotBounds = [(minX - (maxX-minX)/10) (maxX + (maxX-minX)/10) ];
    yPlotBounds = [(minY - (maxY-minY)/10) (maxY + (maxY-minY)/10) ];
end

% label and clean up the plot
axis equal
set(gca,'Ydir','reverse')
title('Distance error')
xlim (xPlotBounds);
ylim (yPlotBounds);

% Create a legend
hSub = subplot(3,3,7);

scatter(nan, nan,2,'filled', ...
    'MarkerFaceAlpha',2/8,'MarkerFaceColor',[0 0 0]);
hold on
scatter(nan, nan,2,'filled', ...
    'MarkerFaceAlpha',2/8,'MarkerFaceColor',[0 0 1]);
plot(nan, nan, '+g', 'MarkerSize', 5);
set(hSub, 'Visible', 'off');
legend({'observed ellipse centers','modeled ellipse centers', 'azimuth 0, elevation 0'},'Location','north', 'Orientation','vertical');


%% Center panel -- shape error
subplot(3,3,[2 5]);

if ~isempty(Xedges)
    % plot the 2D histogram grid
    for xx = 1: length(Xedges)
        if xx==1
            hold on
        end
        plot([Xedges(xx) Xedges(xx)], [Yedges(1) Yedges(end)], '-', 'Color', [0.9 0.9 0.9], 'LineWidth', 0.5 );
    end
    for yy=1: length(Yedges)
        plot([Xedges(1) Xedges(end)], [Yedges(yy) Yedges(yy)], '-', 'Color', [0.9 0.9 0.9], 'LineWidth', 0.5);
    end
end

% Calculate a color for each plot point corresponding to the degree of
% shape error
shapeErrorVec = sceneGeometry.meta.estimateGeometry.search.shapeErrorByEllipse;
shapeErrorVec = shapeErrorVec./sceneGeometry.constraintTolerance;
colorMatrix = zeros(3,size(ellipses,1));
colorMatrix(1,:)=1;
colorMatrix(2,:)= shapeErrorVec;
scatter(ellipses(:,1),ellipses(:,2),[],colorMatrix','o','filled');

% label and clean up the plot
axis equal
set(gca,'Ydir','reverse')
title('Shape error')
xlim (xPlotBounds);
ylim (yPlotBounds);

% Create a legend
hSub = subplot(3,3,8);

scatter(nan, nan,2,'filled', ...
    'MarkerFaceAlpha',6/8,'MarkerFaceColor',[1 0 0]);
hold on
scatter(nan, nan,2,'filled', ...
    'MarkerFaceAlpha',6/8,'MarkerFaceColor',[1 0.5 0]);
scatter(nan, nan,2,'filled', ...
    'MarkerFaceAlpha',6/8,'MarkerFaceColor',[1 1 0]);
set(hSub, 'Visible', 'off');
legend({'0',num2str(sceneGeometry.constraintTolerance/2), ['=> ' num2str(sceneGeometry.constraintTolerance)]},'Location','north', 'Orientation','vertical');


%% Right panel -- area error
subplot(3,3,[3 6]);

if ~isempty(Xedges)
    % plot the 2D histogram grid
    for xx = 1: length(Xedges)
        if xx==1
            hold on
        end
        plot([Xedges(xx) Xedges(xx)], [Yedges(1) Yedges(end)], '-', 'Color', [0.9 0.9 0.9], 'LineWidth', 0.5 );
    end
    for yy=1: length(Yedges)
        plot([Xedges(1) Xedges(end)], [Yedges(yy) Yedges(yy)], '-', 'Color', [0.9 0.9 0.9], 'LineWidth', 0.5);
    end
end

% Calculate a color for each plot point corresponding to the degree of
% shape error
areaErrorVec = sceneGeometry.meta.estimateGeometry.search.areaErrorByEllipse;
areaErrorVec = abs(areaErrorVec)./sceneGeometry.constraintTolerance;
areaErrorVec = min([areaErrorVec ones(size(ellipses,1),1)],[],2);
colorMatrix = zeros(3,size(ellipses,1));
colorMatrix(1,:)=1;
colorMatrix(2,:)= areaErrorVec;
scatter(ellipses(:,1),ellipses(:,2),[],colorMatrix','o','filled');

% label and clean up the plot
axis equal
set(gca,'Ydir','reverse')
title('Area error')
xlim (xPlotBounds);
ylim (yPlotBounds);


% Create a legend
hSub = subplot(3,3,9);

scatter(nan, nan,2,'filled', ...
    'MarkerFaceAlpha',6/8,'MarkerFaceColor',[1 0 0]);
hold on
scatter(nan, nan,2,'filled', ...
    'MarkerFaceAlpha',6/8,'MarkerFaceColor',[1 0.5 0]);
scatter(nan, nan,2,'filled', ...
    'MarkerFaceAlpha',6/8,'MarkerFaceColor',[1 1 0]);
set(hSub, 'Visible', 'off');
legend({'0',num2str(sceneGeometry.constraintTolerance/2), ['=> ' num2str(sceneGeometry.constraintTolerance)]},'Location','north', 'Orientation','vertical');

%% Save the plot
saveas(figHandle,sceneDiagnosticPlotFileName)
close(figHandle)

end % saveSceneDiagnosticPlot


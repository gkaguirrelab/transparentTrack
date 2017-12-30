function initialSceneGeometry = estimateSceneGeometry(pupilFileName, sceneGeometryFileName, varargin)
% Estimate eye radius and eye center given a set of image plane ellipses
%
% Description:
%   This function searches over the set of ellipses in the passed pupil
%   file(s) to estimate the sceneGeometry features that define a
%   perspective projection. The search attempts to minimize the error
%   associated with the prediction of the center of ellipses in the image
%   plane while constraining the shape of these ellipses to match that
%   predicted by the projection.
%
%   We are aware that the center of an ellipse on the image plane is not at
%   the same location as the projection of the center of the pupil on the
%   image plane (see, for example,  Ahn, Sung Joon, H. J. Warnecke, and
%   Rüdiger Kotowski. "Systematic geometric image measurement errors of
%   circular object targets: Mathematical formulation and correction." The
%   Photogrammetric Record 16.93 (1999): 485-502.). The modeling solution
%   implemented here accounts for this property, as we implement a full,
%   numeric forward projection of the pupil circle to the image plane.
%   Similarly, we are aware that there are analytic approximations for the
%   parameters of an ellipse in the image plane given a particular
%   projective circumstance. These analytic approaches may provide an
%   improvement in computation time in comparison to the numeric, forward
%   simulation approach we adopt here, although at some cost in precision.
%
% Components of the projection model:
%   intrinsicCameraMatrix - This matrix has the form:
%
%       [fx  s x0]
%       [0  fy y0]
%       [0   0  1]
%
%   where fx, fy are the focal lengths of the camera in the x and y image
%   dimensions, s is the axis skew, and x0, y0 define the principle offset
%   point. For a camera sensor with square pixels, fx = fy. Ideally, skew
%   should be zero. The principle offset point should be in the center of
%   the sensor. Units are traditionally in pixels. These values can be
%   empirically measured for a camera using a calibration approach
%   (https://www.mathworks.com/help/vision/ref/cameramatrix.html).
%   Therefore, the intrinsicCameraMatrix is locked and not modified in the
%   sceneGeometry search.
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
%   scene coordinate system to be x=0, y=0 at the center of rotation of the
%   eye, and z=0 to be the front surface of the spherical eye ball. The
%   default values and params can be guided by knowledge of the physical
%   arrangement of the subject and recording apparatus.
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
%   eyeRadius -  Scalar which specifies the radius of the eye in mm.
%   Initial value and bounds on the eye radius are taken from:
%
%       Atchison, David A., et al. "Shape of the retinal surface in
%       emmetropia and myopia." Investigative ophthalmology & visual
%       science 46.8 (2005): 2698-2707.
%
%   From Table 1 (mean of axial length and width, and 95% CI)
%           emmetrope:  11.29 (11.07 - 11.51)
%           myope:      11.66 (11.54 - 11.80)
%
%   The default values in the routine are for an emmetrope.
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
% Inputs:
%	pupilFileName         - Full path to a pupilData file, or a cell array
%                           of such paths. If a cell array is provided, the
%                           ellipse data from each pupilData file will be
%                           loaded and concatenated.
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
%  'intrinsicCameraMatrix' - 3x3 matrix
%  'extrinsicRotationMatrix' - 3x3 matrix
%  'eyeRadius'            - Scalar
%  'eyeRadiusLB'          - Scalar
%  'eyeRadiusUB'          - Scalar
%  'extrinsicTranslationVector' - 3x1 matrix
%  'extrinsicTranslationVectorLB' - 3x1 matrix
%  'extrinsicTranslationVectorUB' - 3x1 matrix
%  'constraintTolerance'  - Scalar. Range 0-1.
%  'nBinsPerDimension'    - Scalar. Defines the number of divisions with
%                           which the ellipse centers are binned.
%  'whichEllipseFitField' - Identifies the field in pupilData that contains
%                           the ellipse fit params for which the search
%                           will be conducted.
%
% Outputs
%	sceneGeometry         - A structure with the fields
%       eyeCenter.X - X coordinate of the eye center (i.e. the assumed
%           center of rotation of the pupil) on the scene plane.
%       eyeCenter.Y - Y coordinate of the eye center (i.e. the assumed
%           center of rotation of the pupil) on the scene plane.
%       eyeCenter.Z - the orthogonal distance for the eye center from the
%           scene plane.
%       eyeRadius - radius of the eye in pixels
%       meta - information regarding the analysis, including units.
%


%% input parser
p = inputParser; p.KeepUnmatched = true;

% Required
p.addRequired('pupilFileName',@(x)(iscell(x) | ischar(x)));
p.addRequired('sceneGeometryFileName',@ischar);

% Optional display and I/O params
p.addParameter('verbosity', 'none', @isstr);
p.addParameter('sceneDiagnosticPlotFileName', 'ss',@(x)(isempty(x) | ischar(x)));

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
p.addParameter('intrinsicCameraMatrix',[772.5483 0 320; 0 772.5483 240; 0 0 1],@isnumeric);
p.addParameter('extrinsicRotationMatrix',[1 0 0; 0 -1 0; 0 0 -1],@isnumeric);
p.addParameter('eyeRadius',11.29,@isnumeric);
p.addParameter('eyeRadiusLB',11.07,@isnumeric);
p.addParameter('eyeRadiusUB',11.51,@isnumeric);
p.addParameter('extrinsicTranslationVector',[0; 0; 50],@isnumeric);
p.addParameter('extrinsicTranslationVectorLB',[-10; -10; 45],@isnumeric);
p.addParameter('extrinsicTranslationVectorUB',[10; 10; 65],@isnumeric);
p.addParameter('constraintTolerance',0.02,@isnumeric);
p.addParameter('nBinsPerDimension',7,@isnumeric);
p.addParameter('whichEllipseFitField','initial',@ischar);

% parse
p.parse(pupilFileName, sceneGeometryFileName, varargin{:})


%% Announce we are starting
if strcmp(p.Results.verbosity,'full')
    tic
    fprintf(['Estimating scene geometry from pupil ellipses. Started ' char(datetime('now')) '\n']);
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


%% Identify the ellipses that will guide the sceneGeometry estimation
% load pupil data
if iscell(pupilFileName)
    ellipses = [];
    ellipseFitSEM = [];
    ellipseFitsplitsSD = [];
    for cc = 1:length(pupilFileName)
        load(pupilFileName{cc})
        ellipses = [ellipses;pupilData.(p.Results.whichEllipseFitField).ellipse.values];
        ellipseFitSEM = [ellipseFitSEM; pupilData.(p.Results.whichEllipseFitField).ellipse.RMSE];
    end
else
    load(pupilFileName)
    ellipses = pupilData.(p.Results.whichEllipseFitField).ellipse.values;
    ellipseFitSEM = pupilData.(p.Results.whichEllipseFitField).ellipse.RMSE;
end

% We divide the ellipse centers amongst a set of 2D bins across image
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
[lowestEllipseSEMByBin, idxMinErrorEllipseWithinBin] = arrayfun(@(x) nanmin(ellipseFitSEM(idxByBinPosition{x})), filledBinIdx, 'UniformOutput', false);
errorWeights=cell2mat(lowestEllipseSEMByBin);
errorWeights = 1./errorWeights;
errorWeights=errorWeights./mean(errorWeights);
returnTheMin = @(binContents, x)  binContents(idxMinErrorEllipseWithinBin{x});
ellipseArrayList = cellfun(@(x) returnTheMin(idxByBinPosition{filledBinIdx(x)},x),num2cell(1:1:length(filledBinIdx)));


%% Create the initial sceneGeometry structure and bounds
% sceneGeometry
initialSceneGeometry.eyeRadius = p.Results.eyeRadius;
initialSceneGeometry.intrinsicCameraMatrix = p.Results.intrinsicCameraMatrix;
initialSceneGeometry.extrinsicTranslationVector = p.Results.extrinsicTranslationVector;
initialSceneGeometry.extrinsicRotationMatrix = p.Results.extrinsicRotationMatrix;
initialSceneGeometry.constraintTolerance = p.Results.constraintTolerance;

% Bounds
lb = [p.Results.extrinsicTranslationVectorLB; p.Results.eyeRadiusLB];
ub = [p.Results.extrinsicTranslationVectorUB; p.Results.eyeRadiusUB];


%% Perform the search
% Call out to the local function that performs the serach
sceneGeometry = ...
    performSceneSearch(initialSceneGeometry, ellipses(ellipseArrayList,:), errorWeights, lb, ub);


%% Save the sceneGeometry file
% add a meta field
sceneGeometry.meta = p.Results;
if ~isempty(sceneGeometryFileName)
    save(sceneGeometryFileName,'sceneGeometry');
end


%% Create a sceneGeometry plot
if ~isempty(p.Results.sceneDiagnosticPlotFileName)    
    tmpArray=nan(1,(length(Xedges)-1)*(length(Xedges)-1));
    tmpArray(filledBinIdx)=sceneGeometry.search.errorWeights;
    errorWeightImage = flipud(reshape(tmpArray,length(Xedges)-1,length(Xedges)-1));
    
    tmpArray=nan(1,(length(Xedges)-1)*(length(Xedges)-1));
    tmpArray(filledBinIdx)=sceneGeometry.search.centerDistanceErrorByEllipse;
    centerDistanceErrorByEllipseImage = flipud(reshape(tmpArray,length(Xedges)-1,length(Xedges)-1));
    
    saveSceneDiagnosticPlot(ellipses(ellipseArrayList,:), sceneGeometry.search.errorWeights, Xedges, Yedges, errorWeightImage, centerDistanceErrorByEllipseImage, sceneGeometry, p.Results.sceneDiagnosticPlotFileName)
end

% alert the user that we are done with the routine
if strcmp(p.Results.verbosity,'full')
    toc
    fprintf('\n');
end


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


end % main function



%% LOCAL FUNCTIONS

function sceneGeometry = performSceneSearch(initialSceneGeometry, ellipses, errorWeights, lb, ub)
% Pattern search for best fitting sceneGeometry parameters
%
% Description:
%   The routine searches for parameters of a forward projection that
%   best model the locations of the centers of ellipses found on the image
%   plane, given the constraint that the ellipse shape (and area) must also
%   match the prediction of the forward model. The passed sceneGeometry
%   structure is used as the starting point for the search. The
%   extrinsicTranslationVector and the eyeRadius parameters are optimized,
%   limited by the passed bounds. Across each iteration of the search, a
%   candidate sceneGeometry is assembled from the current values of the
%   parameters. This sceneGeometry is then used in the inverse pupil
%   projection model. The inverse projection searches for an eye azimuth,
%   elevation, and pupil radius that, given the sceneGeometry, best
%   accounts for the parameters of the target ellipse on the image plane.
%   This inverse search attempts to minimize the distance bewteen the
%   centers of the predicted and targeted ellipse on the image plane, while
%   satisfying non-linear constraints upon matching the shape (eccentricity
%   and theta) and area of the ellipses. Only when the sceneGeometry
%   parameters are correctly specified will the inverse pupil projection
%   model be able to simultaneouslty match the center and shape of the
%   ellipse on the image plane.
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

% Set the error norm
norm = 2;

% Extract the initial search point from initialSceneGeometry
x0 = [initialSceneGeometry.extrinsicTranslationVector; initialSceneGeometry.eyeRadius];

% Define search options
options = optimoptions(@patternsearch, ...
    'Display','off',...
    'AccelerateMesh',false,...
    'UseParallel', true, ...
    'FunctionTolerance',0.01);

% Define anonymous functions for the objective and constraint
objectiveFun = @objfun; % the objective function, nested below

% Define nested variables for within the search
centerDistanceErrorByEllipse=[];

[x, fVal] = patternsearch(objectiveFun, x0,[],[],[],[],lb,ub,[],options);
    function fval = objfun(x)
        candidateSceneGeometry = initialSceneGeometry;
        candidateSceneGeometry.extrinsicTranslationVector = x(1:3);
        candidateSceneGeometry.eyeRadius = x(4);
        [~, ~, centerDistanceErrorByEllipse] = ...
            arrayfun(@(x) pupilProjection_inv...
            (...
            ellipses(x,:),...
            candidateSceneGeometry,...
            'constraintTolerance', candidateSceneGeometry.constraintTolerance...
            ),...
            1:1:size(ellipses,1),'UniformOutput',false);
        
        % Now compute objective function as the RMSE of the distance
        % between the taget and modeled ellipses
        centerDistanceErrorByEllipse = cell2mat(centerDistanceErrorByEllipse);
        fval = mean((centerDistanceErrorByEllipse.*errorWeights).^norm)^(1/norm);
    end

% Assemble the sceneGeometry file to return
sceneGeometry.extrinsicTranslationVector = x(1:3);
sceneGeometry.eyeRadius = x(4);
sceneGeometry.intrinsicCameraMatrix = initialSceneGeometry.intrinsicCameraMatrix;
sceneGeometry.extrinsicRotationMatrix = initialSceneGeometry.extrinsicRotationMatrix;
sceneGeometry.constraintTolerance = initialSceneGeometry.constraintTolerance;
sceneGeometry.search.options = options;
sceneGeometry.search.norm = norm;
sceneGeometry.search.initialSceneGeometry = initialSceneGeometry;
sceneGeometry.search.ellipses = ellipses;
sceneGeometry.search.errorWeights = errorWeights;
sceneGeometry.search.lb = lb;
sceneGeometry.search.ub = ub;
sceneGeometry.search.fVal = fVal;
sceneGeometry.search.centerDistanceErrorByEllipse = centerDistanceErrorByEllipse;

end % local search function


function [] = saveSceneDiagnosticPlot(ellipses, errorWeightVec, Xedges, Yedges, errorWeightImage, centerDistanceErrorByEllipseImage, sceneGeometry, sceneDiagnosticPlotFileName)
% Creates and saves a plot that illustrates the sceneGeometry results
%
% Inputs:
%   ellipses              - An n x p array containing the p parameters of
%                           the n ellipses used to derive sceneGeometry
%   errorWeightVec        -
%   Xedges                - The X-dimension edges of the bins used to
%                           divide and select ellipses across the image.
%   Yedges                - The Y-dimension edges of the bins used to
%                           divide and select ellipses across the image.
%   errorWeightImage      - The weight applied to the ellipse at each bin
%                           location in the sceneGeometry search
%   centerDistanceErrorByEllipseImage - The error in the ellipse center
%                           location between the target ellipse and the
%                           best fit of the forward projection model at
%                           each bin location
%   sceneGeometry         - The sceneGeometry structure
%   sceneDiagnosticPlotFileName - The full path (including .pdf suffix)
%                           to the location to save the diagnostic plot
%
% Outputs:
%   none
%

figHandle = figure('visible','off');
subplot(2,2,1)

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
    'constraintTolerance', sceneGeometry.constraintTolerance...
    ),...
    1:1:size(ellipses,1),'UniformOutput',false);
projectedEllipses=vertcat(projectedEllipses{:});

% plot the projected ellipse centers
scatter(projectedEllipses(:,1),projectedEllipses(:,2),'o','filled', ...
    'MarkerFaceAlpha',2/8,'MarkerFaceColor',[0 0 1]);

% connect the centers with lines
for ii=1:size(ellipses,1)
    lineAlpha = errorWeightVec(ii)/max(errorWeightVec);
    lineWeight = 0.5 + (errorWeightVec(ii)/max(errorWeightVec));
    ph=plot([projectedEllipses(ii,1) ellipses(ii,1)], ...
        [projectedEllipses(ii,2) ellipses(ii,2)], ...
        '-','Color',[1 0 0],'LineWidth', lineWeight);
    ph.Color(4) = lineAlpha;
end

% plot the estimated center of rotation of the eye
centerOfRotationEllipse = pupilProjection_fwd([0 0 2], sceneGeometry);
plot(centerOfRotationEllipse(1),centerOfRotationEllipse(2), '+g', 'MarkerSize', 5);

% label and clean up the plot
xlim ([Xedges(1)-binSpaceX Xedges(end)+binSpaceX]);
ylim ([Yedges(1)-binSpaceY Yedges(end)+binSpaceY]);
axis equal
set(gca,'Ydir','reverse')
title('Ellipse centers')

% Create a legend
hSub = subplot(2,2,2);
scatter(nan, nan,2,'filled', ...
    'MarkerFaceAlpha',2/8,'MarkerFaceColor',[0 0 0]);
hold on
scatter(nan, nan,2,'filled', ...
    'MarkerFaceAlpha',2/8,'MarkerFaceColor',[0 0 1]);
plot(nan, nan, '+g', 'MarkerSize', 5);
set(hSub, 'Visible', 'off');
legend({'observed ellipse centers','modeled ellipse centers', 'azimuth 0, elevation 0'},'Location','southwestoutside');

% Plot the ellipse counts and error values by bin
subplot(2,2,3)
nanAwareImagePlot(errorWeightImage, Xedges, Yedges, 'Error weights')

% Plot the centerDistanceErrorByEllipse in each bin
subplot(2,2,4)
nanAwareImagePlot(centerDistanceErrorByEllipseImage, Xedges, Yedges, 'Ellipse center distance error')

% Save the plot
saveas(figHandle,sceneDiagnosticPlotFileName);
close(figHandle)

end % saveSceneDiagnosticPlot


function nanAwareImagePlot(image, Xedges, Yedges, titleString)
% A replacement for imagesc that gracefully handles nans
%
[nr,nc] = size(image);
pcolor([flipud(image) nan(nr,1); nan(1,nc+1)]);
caxis([0 max(max(image))]);
shading flat;
axis equal

% Set the axis backgroud to dark gray
set(gcf,'Color',[1 1 1]); set(gca,'Color',[.75 .75 .75]); set(gcf,'InvertHardCopy','off');
set(gca,'Ydir','reverse')
colorbar;
title(titleString);
xticks(1:1:size(image,1)+1);
xticklabels(round(Xedges));
xtickangle(90);
yticks(1:1:size(image,2)+1);
yticklabels(round(Yedges));
xlim([1 size(image,1)+1]);
ylim([1 size(image,2)+1]);
end

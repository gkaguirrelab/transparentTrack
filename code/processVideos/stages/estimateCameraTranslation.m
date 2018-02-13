function sceneGeometry = estimateCameraTranslation(pupilFileName, sceneGeometryFileName, varargin)
% Estimate camera translation given a set of image plane ellipses
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
%   forward projection of the pupil circle to the image plane.
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
%  'translationLB/UB'     - 3x1 vector. Hard upper and lower bounds. Should
%                           reflect the physical limits of the measurement.
%  'translationLBp/UBp'   - 3x1 vector. Plausible upper and lower bounds.
%                           Where you think the translation vector solution
%                           is likely to be.
%  'eyePoseLB/UB'         - 1x4 vector. Upper / lower bounds on the eyePose
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
%                           which identify the ellipses to be used for the
%                           estimation of scene geometry. If left empty,
%                           a list of ellipses will be generated.
%  'nBinsPerDimension'    - Scalar. Defines the number of divisions with
%                           which the ellipse centers are binned.
%  'useRayTracing'        - Logical; default false. Using ray tracing in
%                           the camera translation search improves accuracy
%                           slightly, but increases search time by about
%                           25x.
%  'nBADSsearches'        - Scalar. We perform the search for camera 
%                           translation from a randomly selected starting
%                           point within the plausible bounds. This
%                           parameter sets how many random starting points
%                           to try; the best result is retained. Each
%                           search is run on a separate worker if the
%                           parpool is available.
%                           
%
% Outputs
%	sceneGeometry         - A structure that contains the components of the
%                           projection model.
%
% Examples:
%{
    %% Recover a veridical camera translation
    % Create a veridical sceneGeometry with some arbitrary translation
    veridicalSceneGeometry = createSceneGeometry();
    veridicalSceneGeometry.extrinsicTranslationVector = [-1.2; 0.9; 108];
    % Assemble the ray tracing functions
    rayTraceFuncs = assembleRayTraceFuncs( veridicalSceneGeometry );
    % Create a set of ellipses using the veridical geometry and 
    % randomly varying pupil radii.
    ellipseIdx=1;
    for azi=-15:15:15
    	for ele=-15:15:15
            eyePose=[azi, ele, 0, 2+(randn()./5)];
            pupilData.initial.ellipses.values(ellipseIdx,:) = pupilProjection_fwd(eyePose, veridicalSceneGeometry, rayTraceFuncs);
            pupilData.initial.ellipses.RMSE(ellipseIdx,:) = 1;
            ellipseIdx=ellipseIdx+1;
        end
    end
    % Estimate the scene Geometry using the ellipses
    estimatedSceneGeometry = estimateCameraTranslation(pupilData,'','useParallel',true,'verbosity','full','ellipseArrayList',1:1:ellipseIdx-1,'nBADSsearches',2);
    % Report how well we did
    fprintf('Error in the recovered camera translation vector (x, y, depth] in mm: \n');
    veridicalSceneGeometry.extrinsicTranslationVector - estimatedSceneGeometry.extrinsicTranslationVector
%}

%% input parser
p = inputParser; p.KeepUnmatched = true;

% Required
p.addRequired('pupilFileName',@(x)(isempty(x) | isstruct(x) | iscell(x) | ischar(x)));
p.addRequired('sceneGeometryFileName',@ischar);

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
p.addParameter('translationLB',[-20; -20; 90],@isnumeric);
p.addParameter('translationUB',[20; 20; 200],@isnumeric);
p.addParameter('translationLBp',[-10; -10; 100],@isnumeric);
p.addParameter('translationUBp',[5; 5; 160],@isnumeric);
p.addParameter('eyePoseLB',[-35,-25,0,0.25],@(x)(isempty(x) | isnumeric(x)));
p.addParameter('eyePoseUB',[35,25,0,4],@(x)(isempty(x) | isnumeric(x)));
p.addParameter('fitLabel','initial',@ischar);
p.addParameter('ellipseArrayList',[],@(x)(isempty(x) | isnumeric(x)));
p.addParameter('nBinsPerDimension',4,@isnumeric);
p.addParameter('useRayTracing',false,@islogical);
p.addParameter('nBADSsearches',10,@isnumeric);

% parse
p.parse(pupilFileName, sceneGeometryFileName, varargin{:})


%% Announce we are starting
if strcmp(p.Results.verbosity,'full')
    tic
    fprintf(['Estimating camera translation from pupil ellipses. Started ' char(datetime('now')) '\n']);
end

%% Create initial sceneGeometry structure and ray tracing functions
initialSceneGeometry = createSceneGeometry(varargin{:});

% Assemble the ray tracing functions
if p.Results.useRayTracing
    if strcmp(p.Results.verbosity,'full')
        fprintf('Assembling ray tracing functions.\n');
    end
    [rayTraceFuncs] = assembleRayTraceFuncs( initialSceneGeometry );
else
    rayTraceFuncs = [];
end

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
        ellipses = [ellipses; pupilData.(p.Results.fitLabel).ellipses.values];
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
    
    % First we divide the ellipse centers amongst a set of 2D bins across
    % image space.
    [ellipseCenterCounts,Xedges,Yedges,binXidx,binYidx] = ...
        histcounts2(ellipses(:,1),ellipses(:,2),p.Results.nBinsPerDimension);
    
    % Anonymous functions for row and column identity given array position
    rowIdx = @(b) fix( (b-1) ./ (size(ellipseCenterCounts,2)) ) +1;
    colIdx = @(b) 1+mod(b-1,size(ellipseCenterCounts,2));
    
    % Create a cell array of index positions corresponding to each of the
    % 2D bins
    idxByBinPosition = ...
        arrayfun(@(b) find( (binXidx==rowIdx(b)) .* (binYidx==colIdx(b)) ),1:1:numel(ellipseCenterCounts),'UniformOutput',false);
    
    % Identify which bins are not empty
    filledBinIdx = find(~cellfun(@isempty, idxByBinPosition));
    
    % Identify the ellipse in each bin with the lowest fit SEM
    [~, idxMinErrorEllipseWithinBin] = arrayfun(@(x) nanmin(ellipseFitSEM(idxByBinPosition{x})), filledBinIdx, 'UniformOutput', false);
    returnTheMin = @(binContents, x)  binContents(idxMinErrorEllipseWithinBin{x});
    ellipseArrayList = cellfun(@(x) returnTheMin(idxByBinPosition{filledBinIdx(x)},x),num2cell(1:1:length(filledBinIdx)));
end


%% Generate the errorWeights
errorWeights = ellipseFitSEM(ellipseArrayList);
errorWeights = 1./errorWeights;
errorWeights = errorWeights./mean(errorWeights);


%% Perform the search
if strcmp(p.Results.verbosity,'full')
    fprintf(['Searching over camera translations.\n']);
    fprintf('| 0                      50                   100%% |\n');
    fprintf('.\n');
end

searchResults = {};
parfor (ss = 1:p.Results.nBADSsearches,nWorkers)
    
    searchResults{ss} = ...
        performSceneSearch(initialSceneGeometry, rayTraceFuncs, ...
        ellipses(ellipseArrayList,:), ...
        errorWeights, ...
        p.Results.translationLB, ...
        p.Results.translationUB, ...
        p.Results.translationLBp, ...
        p.Results.translationUBp, ...
        p.Results.eyePoseLB, ...
        p.Results.eyePoseUB);
    
    % update progress
    if strcmp(p.Results.verbosity,'full')
        for pp=1:floor(50/p.Results.nBADSsearches)
            fprintf('\b.\n');
        end
    end
    
end
if strcmp(p.Results.verbosity,'full')
    fprintf('\n');
end

% Keep the best result
allFvals = cellfun(@(x) x.meta.estimateCameraTranslation.search.fVal,searchResults);
allTranslationVecs = cellfun(@(x) x.extrinsicTranslationVector,searchResults,'UniformOutput',false);
[~,idx]=min(allFvals);
sceneGeometry = searchResults{idx};
    
% add additional search and meta field info to sceneGeometry
sceneGeometry.meta.estimateCameraTranslation.parameters = p.Results;
sceneGeometry.meta.estimateCameraTranslation.search.ellipseArrayList = ellipseArrayList';
sceneGeometry.meta.estimateCameraTranslation.search.allFvals = allFvals;
sceneGeometry.meta.estimateCameraTranslation.search.allTranslationVecs = allTranslationVecs;


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
        p.Results.eyePoseLB, ...
        p.Results.eyePoseUB, ...
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

function sceneGeometry = performSceneSearch(initialSceneGeometry, rayTraceFuncs, ellipses, errorWeights, LB, UB, LBp, UBp, eyePoseLB, eyePoseUB)
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
%   The search is performed using Bayesian Adaptive Direct Search (bads),
%   as we find that it performs better than (e.g.) patternsearch. BADS only
%   accepts row vectors, so there is much transposing ahead.
%

% Set the error form
errorForm = 'RMSE';

% Pick a random x0 from within the plausible bounds
x0 = LBp + (UBp-LBp).*rand(numel(LBp),1);

% Define search options
options = bads('defaults');          % Get a default OPTIONS struct
options.Display = 'off';             % Silence display output
options.UncertaintyHandling = 0;     % The objective is deterministic

% Silence the mesh overflow warning from BADS
warningState = warning;
warning('off','bads:meshOverflow');

% Define nested variables for within the search
centerDistanceErrorByEllipse=zeros(size(ellipses,1),1);
shapeErrorByEllipse=zeros(size(ellipses,1),1);
areaErrorByEllipse=zeros(size(ellipses,1),1);

% Perform the seach using bads
[x, fVal] = bads(@objfun,x0',LB',UB',LBp',UBp',[],options);
% Nested function computes the objective
    function fval = objfun(x)
        % Assemble a candidate sceneGeometry structure
        candidateSceneGeometry = initialSceneGeometry;
        candidateSceneGeometry.extrinsicTranslationVector = x';
        % For each ellipse, perform the inverse projection from the ellipse
        % on the image plane to eyePose. We retain the errors from the
        % inverse projection and use these to assemble the objective
        % function. We parallelize the computation across ellipses.
        for ii = 1:size(ellipses,1)
            exitFlag = [];
            eyePose = [];
            [eyePose, ~, centerDistanceErrorByEllipse(ii), shapeErrorByEllipse(ii), areaErrorByEllipse(ii), exitFlag] = ...
                pupilProjection_inv(...
                ellipses(ii,:),...
                candidateSceneGeometry, rayTraceFuncs, ...
                'eyePoseLB',eyePoseLB,...
                'eyePoseUB',eyePoseUB...
                );
            % if the exitFlag indicates a possible local minimum, repeat
            % the search and initialize with the returned eyePose
            if exitFlag == 2
                x0tmp = eyePose + [1e-3 1e-3 0 1-3];
                [~, ~, centerDistanceErrorByEllipse(ii), shapeErrorByEllipse(ii), areaErrorByEllipse(ii)] = ...
                    pupilProjection_inv(...
                    ellipses(ii,:),...
                    candidateSceneGeometry, rayTraceFuncs, ...
                    'eyePoseLB',eyePoseLB,...
                    'eyePoseUB',eyePoseUB,...
                    'x0',x0tmp...
                    );
            end
        end
        % Now compute objective function as the RMSE of the distance
        % between the taget and modeled ellipses
        switch errorForm
            case 'SSE'
                fval=sum((centerDistanceErrorByEllipse.*(shapeErrorByEllipse.*100+1).*(areaErrorByEllipse.*100+1).*errorWeights).^2);
                % We have to keep the fval non-infinite to keep bads happy
                fval=min([fval realmax]);
            case 'RMSE'
                fval = mean((centerDistanceErrorByEllipse.*(shapeErrorByEllipse.*100+1).*(areaErrorByEllipse.*100+1).*errorWeights).^2).^(1/2);
                % We have to keep the fval non-infinite to keep bads happy
                fval=min([fval realmax]);
            otherwise
                error('I do not recognize that error form');
        end
        
    end

% Restore the warning state
warning(warningState);

% Assemble the sceneGeometry file to return
sceneGeometry.radialDistortionVector = initialSceneGeometry.radialDistortionVector;
sceneGeometry.intrinsicCameraMatrix = initialSceneGeometry.intrinsicCameraMatrix;
sceneGeometry.extrinsicTranslationVector = x';
sceneGeometry.extrinsicRotationMatrix = initialSceneGeometry.extrinsicRotationMatrix;
sceneGeometry.primaryPosition = initialSceneGeometry.primaryPosition;
sceneGeometry.constraintTolerance = initialSceneGeometry.constraintTolerance;
sceneGeometry.eye = initialSceneGeometry.eye;
sceneGeometry.meta.estimateCameraTranslation.search.options = options;
sceneGeometry.meta.estimateCameraTranslation.search.errorForm = errorForm;
sceneGeometry.meta.estimateCameraTranslation.search.initialSceneGeometry = initialSceneGeometry;
sceneGeometry.meta.estimateCameraTranslation.search.ellipses = ellipses;
sceneGeometry.meta.estimateCameraTranslation.search.errorWeights = errorWeights;
sceneGeometry.meta.estimateCameraTranslation.search.x0 = x0;
sceneGeometry.meta.estimateCameraTranslation.search.LB = LB;
sceneGeometry.meta.estimateCameraTranslation.search.UB = UB;
sceneGeometry.meta.estimateCameraTranslation.search.LBp = LBp;
sceneGeometry.meta.estimateCameraTranslation.search.UBp = UBp;
sceneGeometry.meta.estimateCameraTranslation.search.eyePoseLB = eyePoseLB;
sceneGeometry.meta.estimateCameraTranslation.search.eyePoseUB = eyePoseUB;
sceneGeometry.meta.estimateCameraTranslation.search.fVal = fVal;
sceneGeometry.meta.estimateCameraTranslation.search.centerDistanceErrorByEllipse = centerDistanceErrorByEllipse;
sceneGeometry.meta.estimateCameraTranslation.search.shapeErrorByEllipse = shapeErrorByEllipse;
sceneGeometry.meta.estimateCameraTranslation.search.areaErrorByEllipse = areaErrorByEllipse;

end % local search function


function [] = saveSceneDiagnosticPlot(ellipses, Xedges, Yedges, eyePoseLB, eyePoseUB, sceneGeometry, rayTraceFuncs, sceneDiagnosticPlotFileName)
% Creates and saves a plot that illustrates the sceneGeometry results
%
% Inputs:
%   ellipses              - An n x p array containing the p parameters of
%                           the n ellipses used to derive sceneGeometry
%   Xedges                - The X-dimension edges of the bins used to
%                           divide and select ellipses across the image.
%   Yedges                - The Y-dimension edges of the bins used to
%                           divide and select ellipses across the image.
%   eyePoseLB, eyePoseUB  - Bounds for the eye pose to be passed to
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
    'Renderer','painters'...
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
    'eyePoseLB',eyePoseLB,'eyePoseUB',eyePoseUB),...
    1:1:size(ellipses,1),'UniformOutput',false);
projectedEllipses=vertcat(projectedEllipses{:});

% plot the projected ellipse centers
scatter(projectedEllipses(:,1),projectedEllipses(:,2),'o','filled', ...
    'MarkerFaceAlpha',2/8,'MarkerFaceColor',[0 0 1]);

% connect the centers with lines
errorWeightVec=sceneGeometry.meta.estimateCameraTranslation.search.errorWeights;
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
shapeErrorVec = sceneGeometry.meta.estimateCameraTranslation.search.shapeErrorByEllipse;
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
areaErrorVec = sceneGeometry.meta.estimateCameraTranslation.search.areaErrorByEllipse;
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


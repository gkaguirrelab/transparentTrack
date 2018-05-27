function sceneGeometry = estimateSceneParams(pupilFileName, sceneGeometryFileName, varargin)
% Estimate camera translation and eye rotation given image plane ellipses
%
% Syntax:
%  sceneGeometry = estimateSceneParams(pupilFileName, sceneGeometryFileName)
%
% Description:
%   This function searches over a set of ellipses from the passed pupil
%   file(s) to estimate the extrinsic camera rotation and translation
%   vector and scaling values for the azimuthal and elevational eye
%   rotation centers. The search attempts to minimize the error associated
%   with the prediction of the shape of ellipses in the image plane while
%   minimizing the error in prediction of the center of those ellipses in
%   the image plane.
%
%   The search is conducted over 6 parameters, corresponding to camera
%   rotation about the Z (depth) axis, three parameters of camera
%   translation (horizontal, vertical, depth), a parameter for joint
%   scaling of the centers of rotation of the eye (azimuthal and
%   elevational rotations), and then a parameter for differential scaling
%   of the eye rotation centers. For this last parameter, a value > 1
%   increases the azimuthal rotation center values and decreases the
%   elevational.
%
% Inputs:
%	pupilFileName         - Full path to a pupilData file, a cell array
%                           of such paths, or a pupilData structure itself.
%                           If a single path, the pupilData file is loaded.
%                           If a cell array, the ellipse data from each
%                           pupilData file is loaded and concatenated.
%   sceneGeometryFileName - Full path to the file in which the
%                           sceneGeometry data should be saved
%   vitualImageFuncDir    - Full path to the directory that should be
%                           created to hold the ray tracing function for
%                           the optical system in sceneGeometry.
%
% Optional key/value pairs (display and I/O):
%  'verbose'              - Logical. Default false.
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
%  'sceneParamsLB/UB'     - 6x1 vector. Hard upper and lower bounds. Should
%                           reflect the physical limits of the measurement.
%  'sceneParamsLBp/UBp'   - 6x1 vector. Plausible upper and lower bounds.
%                           Where you think the translation vector solution
%                           is likely to be.
%  'eyePoseLB/UB'         - 1x4 vector. Upper / lower bounds on the eyePose
%                           [azimuth, elevation, torsion, pupil radius].
%                           The torsion value is unusued and is bounded to
%                           zero. Biological limits in eye rotation and
%                           pupil size would suggest boundaries of [±35,
%                           ±25, 0, 0.25-4]. Note, however, that these
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
%  'badFrameErrorThreshold' - Frames with RMSE fitting error above this
%                           threshold will not be selected to guide the
%                           scene parameter search.
%  'nBADSsearches'        - Scalar. We perform the search over scene params
%                           from a randomly selected starting point within
%                           the plausible bounds. This parameter sets how
%                           many random starting points to try; the best
%                           result is retained. Each search is run on a
%                           separate worker if the parpool is available.
%
%
% Outputs
%	sceneGeometry         - A structure that contains the components of the
%                           projection model.
%
% Examples:
%{
    %% Recover a veridical camera translation
    % Create a default sceneGeometry
    veridicalSceneGeometry = createSceneGeometry();
    % Set some arbitrary extrinsic camera translation
    veridicalSceneGeometry.cameraPosition.translation = [-1.2; 0.9; 108];
    % Create a set of ellipses using the veridical geometry and
    % randomly varying pupil radii.
    ellipseIdx=1;
    for azi=-15:15:15
    	for ele=-15:15:15
            eyePose=[azi, ele, 0, 2+(randn()./5)];
            pupilData.initial.ellipses.values(ellipseIdx,:) = pupilProjection_fwd(eyePose, veridicalSceneGeometry);
            pupilData.initial.ellipses.RMSE(ellipseIdx,:) = 1;
            ellipseIdx=ellipseIdx+1;
        end
    end
    % Estimate the scene Geometry using the ellipses
    estimatedSceneGeometry = estimateSceneParams(pupilData,'','useParallel',true,'verbose',true,'ellipseArrayList',1:1:ellipseIdx-1);
    % Report how well we did
    fprintf('Error in the recovered camera translation vector (x, y, depth] in mm: \n');
    veridicalSceneGeometry.cameraPosition.translation - estimatedSceneGeometry.cameraPosition.translation
%}

%% input parser
p = inputParser; p.KeepUnmatched = true;

% Required
p.addRequired('pupilFileName',@(x)(isstruct(x) || iscell(x) || ischar(x)));
p.addRequired('sceneGeometryFileName',@ischar);

% Optional display and I/O params
p.addParameter('verbose',false,@islogical);
p.addParameter('sceneDiagnosticPlotFileName', '', @(x)(isempty(x) || ischar(x)));

% Optional flow control params
p.addParameter('useParallel',false,@islogical);
p.addParameter('nWorkers',[],@(x)(isempty(x) || isnumeric(x)));

% Optional environment parameters
p.addParameter('tbSnapshot',[],@(x)(isempty(x) || isstruct(x)));
p.addParameter('timestamp',char(datetime('now')),@ischar);
p.addParameter('username',char(java.lang.System.getProperty('user.name')),@ischar);
p.addParameter('hostname',char(java.net.InetAddress.getLocalHost.getHostName),@ischar);

% Optional analysis params
p.addParameter('sceneParamsLB',[-30; -20; -20; 90; 0.75; .9],@isnumeric);
p.addParameter('sceneParamsUB',[30; 20; 20; 200; 1.25; 1.1],@isnumeric);
p.addParameter('sceneParamsLBp',[-15; -5; -5; 100; 0.85; 0.95],@isnumeric);
p.addParameter('sceneParamsUBp',[15; 5; 5; 160; 1.15; 1.05],@isnumeric);
p.addParameter('eyePoseLB',[-89,-89,0,0.1],@isnumeric);
p.addParameter('eyePoseUB',[89,89,0,4],@isnumeric);
p.addParameter('fitLabel','initial',@ischar);
p.addParameter('ellipseArrayList',[],@(x)(isempty(x) | isnumeric(x)));
p.addParameter('nBinsPerDimension',4,@isnumeric);
p.addParameter('badFrameErrorThreshold',2, @isnumeric);
p.addParameter('nBADSsearches',10,@isnumeric);

% parse
p.parse(pupilFileName, sceneGeometryFileName, varargin{:})


%% Announce we are starting
if p.Results.verbose
    tic
    fprintf(['Estimating camera position and eye rotation from pupil ellipses. Started ' char(datetime('now')) '\n']);
end

%% Create initial sceneGeometry structure
initialSceneGeometry = createSceneGeometry(varargin{:});


%% Set up the parallel pool
if p.Results.useParallel
    nWorkers = startParpool( p.Results.nWorkers, p.Results.verbose );
else
    nWorkers=0;
end


%% Load pupil data
if iscell(pupilFileName)
    ellipses = [];
    ellipseFitRMSE = [];
    for cc = 1:length(pupilFileName)
        load(pupilFileName{cc})
        ellipses = [ellipses; pupilData.(p.Results.fitLabel).ellipses.values];
        ellipseFitRMSE = [ellipseFitRMSE; pupilData.(p.Results.fitLabel).ellipses.RMSE];
    end
end
if ischar(pupilFileName)
    load(pupilFileName)
    ellipses = pupilData.(p.Results.fitLabel).ellipses.values;
    ellipseFitRMSE = pupilData.(p.Results.fitLabel).ellipses.RMSE;
end
if isstruct(pupilFileName)
    pupilData = pupilFileName;
    ellipses = pupilData.(p.Results.fitLabel).ellipses.values;
    ellipseFitRMSE = pupilData.(p.Results.fitLabel).ellipses.RMSE;
end


%% Identify the ellipses that will guide the sceneGeometry estimation
% If not supplied, we will generate a list of ellipses to use for the
% estimation.
if ~isempty(p.Results.ellipseArrayList)
    ellipseArrayList = p.Results.ellipseArrayList;
    Xedges = [];
    Yedges = [];
else
    if p.Results.verbose
        fprintf('Selecting ellipses to guide the search.\n');
    end
    
    % Identify the ellipses with RMSE fits that are below the "bad"
    % threshold
    goodFitIdx = find(ellipseFitRMSE < p.Results.badFrameErrorThreshold);
    if isempty(goodFitIdx)
        error('No initial ellipse fits are good enough to guide the search; try adjusting badFrameErrorThreshold');
    end
    
    % First we divide the ellipse centers amongst a set of 2D bins across
    % image space.
    [ellipseCenterCounts,Xedges,Yedges,binXidx,binYidx] = ...
        histcounts2(ellipses(goodFitIdx,1),ellipses(goodFitIdx,2),p.Results.nBinsPerDimension);
    
    % Anonymous functions for row and column identity given array position
    rowIdx = @(b) fix( (b-1) ./ (size(ellipseCenterCounts,2)) ) +1;
    colIdx = @(b) 1+mod(b-1,size(ellipseCenterCounts,2));
    
    % Create a cell array of index positions corresponding to each of the
    % 2D bins
    idxByBinPosition = ...
        arrayfun(@(b) find( (binXidx==rowIdx(b)) .* (binYidx==colIdx(b)) ),1:1:numel(ellipseCenterCounts),'UniformOutput',false);
    
    % Identify which bins are not empty
    filledBinIdx = find(~cellfun(@isempty, idxByBinPosition));
    
    % Identify the ellipse in each bin with the lowest fit RMSE
    [~, idxMinErrorEllipseWithinBin] = arrayfun(@(x) nanmin(ellipseFitRMSE(goodFitIdx(idxByBinPosition{x}))), filledBinIdx, 'UniformOutput', false);
    returnTheMin = @(binContents, x)  binContents(idxMinErrorEllipseWithinBin{x});
    ellipseArrayList = cellfun(@(x) returnTheMin(goodFitIdx(idxByBinPosition{filledBinIdx(x)}),x),num2cell(1:1:length(filledBinIdx)));
end


%% Perform the search
% Inform the user
if p.Results.verbose
    fprintf(['Searching over scene geometry parameters.\n']);
    fprintf('| 0                      50                   100%% |\n');
    fprintf('.\n');
end

% Loop over the requested number of BADS searches
searchResults = {};
parfor (ss = 1:p.Results.nBADSsearches,nWorkers)
    
    searchResults{ss} = ...
        performSceneSearch(initialSceneGeometry, ...
        ellipses(ellipseArrayList,:), ...
        p.Results.sceneParamsLB, ...
        p.Results.sceneParamsUB, ...
        p.Results.sceneParamsLBp, ...
        p.Results.sceneParamsUBp, ...
        p.Results.eyePoseLB, ...
        p.Results.eyePoseUB, ...
        p.Results.badFrameErrorThreshold);
    
    % update progress
    if p.Results.verbose
        for pp=1:floor(50/p.Results.nBADSsearches(1))
            fprintf('\b.\n');
        end
    end
    
end
if p.Results.verbose
    fprintf('\n');
end

% Find the weighted mean and SD of the solution parameters
allFvals = cellfun(@(x) x.meta.estimateSceneParams.search.fVal,searchResults);
allSceneParamResults = cellfun(@(thisSceneGeometry) thisSceneGeometry.meta.estimateSceneParams.search.x,searchResults,'UniformOutput',false);
for dim = 1:length(p.Results.sceneParamsLB)
    vals = cellfun(@(x) x(dim), allSceneParamResults);
    sceneParamResultsMean(dim)=mean(vals.*(1./allFvals))/mean(1./allFvals);
    sceneParamResultsSD(dim)=std(vals,1./allFvals);
end
sceneParamResultsMean=sceneParamResultsMean';
sceneParamResultsSD=sceneParamResultsSD';

% Perform the search using the mean parameters as absolute bounds to obtain
% the error values
sceneGeometry = ...
    performSceneSearch(initialSceneGeometry, ...
    ellipses(ellipseArrayList,:), ...
    sceneParamResultsMean, ...
    sceneParamResultsMean, ...
    sceneParamResultsMean, ...
    sceneParamResultsMean, ...
    p.Results.eyePoseLB, ...
    p.Results.eyePoseUB);

% Add additional search and meta field info to sceneGeometry
tmpHold=sceneGeometry.meta.estimateSceneParams.search;
sceneGeometry.meta.estimateSceneParams = p.Results;
sceneGeometry.meta.estimateSceneParams.search = tmpHold;
sceneGeometry.meta.estimateSceneParams.search.ellipseArrayList = ellipseArrayList';
sceneGeometry.meta.estimateSceneParams.search.ellipseRMSE = ellipseFitRMSE(ellipseArrayList);
sceneGeometry.meta.estimateSceneParams.search.allFvals = allFvals;
sceneGeometry.meta.estimateSceneParams.search.allSceneParamResults = allSceneParamResults;
sceneGeometry.meta.estimateSceneParams.search.sceneParamResultsMean = sceneParamResultsMean;
sceneGeometry.meta.estimateSceneParams.search.sceneParamResultsSD = sceneParamResultsSD;


%% Save the sceneGeometry file
if ~isempty(sceneGeometryFileName)
    save(sceneGeometryFileName,'sceneGeometry');
end


%% Create a sceneGeometry plot
if ~isempty(p.Results.sceneDiagnosticPlotFileName)
    if p.Results.verbose
        fprintf('Creating a sceneGeometry diagnostic plot.\n');
    end
    saveSceneDiagnosticPlot(...
        ellipses(ellipseArrayList,:),...
        Xedges, Yedges,...
        p.Results.eyePoseLB, ...
        p.Results.eyePoseUB, ...
        sceneGeometry,...
        p.Results.sceneDiagnosticPlotFileName)
end


%% alert the user that we are done with the routine
if p.Results.verbose
    toc
    fprintf('\n');
end


end % main function



%% LOCAL FUNCTIONS

function sceneGeometry = performSceneSearch(initialSceneGeometry, ellipses, LB, UB, LBp, UBp, eyePoseLB, eyePoseUB, badFrameErrorThreshold)
% Pattern search for best fitting sceneGeometry parameters
%
% Description:
%   The routine searches for parameters of camera translation and torsion,
%   and rotation centers of the eye, that best model the shapes (and areas)
%   of ellipses found on the image plane, while minimizing the distance
%   between the modeled and observed ellipse centers. The passed
%   sceneGeometry structure is used as the starting point for the search.
%   Across each iteration of the search, a candidate sceneGeometry is
%   assembled from the current values of the parameters. This sceneGeometry
%   is then used in the inverse pupil projection model. The inverse
%   projection searches for an eye azimuth, elevation, and pupil radius
%   that, given the sceneGeometry, best accounts for the parameters of the
%   target ellipse on the image plane. This inverse search attempts to
%   minimize the distance bewteen the centers of the predicted and targeted
%   ellipse on the image plane, while satisfying non-linear constraints
%   upon matching the shape (eccentricity and theta) and area of the
%   ellipses. Only when the translation vector is correctly specified will
%   the inverse pupil projection model be able to simultaneouslty match the
%   center and shape of the ellipse on the image plane.
%
%   The iterative search across sceneGeometry parameters attempts to
%   minimize the L2 norm of the shape and area errors between the targeted
%   and modeled centers of the ellipses.
%
%   The search is performed using Bayesian Adaptive Direct Search (bads),
%   as we find that it performs better than (e.g.) patternsearch. BADS only
%   accepts row vectors, so there is much transposing ahead.
%

% Pick a random x0 from within the plausible bounds
x0 = LBp + (UBp-LBp).*rand(numel(LBp),1);

% Define search options
options = bads('defaults');          % Get a default OPTIONS struct
options.Display = 'off';             % Silence display output
options.UncertaintyHandling = 0;     % The objective is deterministic

% Silence the mesh overflow warning from BADS
warningState = warning;
warning('off','bads:meshOverflow');

% Silence some errors that can arise during the forward projection
warning('off','rayTraceEllipsoids:criticalAngle');
warning('off','pupilProjection_fwd:ellipseFitFailed');

% Define nested variables for within the search
centerDistanceErrorByEllipse=zeros(size(ellipses,1),1);
shapeErrorByEllipse=zeros(size(ellipses,1),1);
areaErrorByEllipse=zeros(size(ellipses,1),1);
recoveredEyePoses =zeros(size(ellipses,1),4);

% Detect if we have pinned the parameters, in which case just evaluate the
% objective function
if all(x0==LB) && all(x0==UB)
    x=x0';
    fVal = objfun(x);
else
    % Perform the seach using bads
    [x, fVal] = bads(@objfun,x0',LB',UB',LBp',UBp',[],options);
end
% Nested function computes the objective
    function fval = objfun(x)
        % Assemble a candidate sceneGeometry structure
        candidateSceneGeometry = initialSceneGeometry;
        % Store the camera torsion
        candidateSceneGeometry.cameraPosition.torsion = x(1);
        % Store the extrinsic camera translation vector
        candidateSceneGeometry.cameraPosition.translation = x(2:4)';
        % Scale the rotation center values by the joint and differential
        % parameters
        candidateSceneGeometry.eye.rotationCenters.azi = candidateSceneGeometry.eye.rotationCenters.azi .* x(5) .* x(6);
        candidateSceneGeometry.eye.rotationCenters.ele = candidateSceneGeometry.eye.rotationCenters.ele .* x(5) ./ x(6);
        % For each ellipse, perform the inverse projection from the ellipse
        % on the image plane to eyePose. We retain the errors from the
        % inverse projection and use these to assemble the objective
        % function.
        for ii = 1:size(ellipses,1)
            [recoveredEyePoses(ii,:), ~, centerDistanceErrorByEllipse(ii), shapeErrorByEllipse(ii), areaErrorByEllipse(ii)] = ...
                pupilProjection_inv(...
                    ellipses(ii,:),...
                    candidateSceneGeometry, ...
                    'eyePoseLB',eyePoseLB,...
                    'eyePoseUB',eyePoseUB,...
                    'repeatSearchThresh',badFrameErrorThreshold);
        end
        % Now compute objective function as the RMSE of the distance
        % between the taget and modeled ellipses in shape and area
        fval = sqrt(mean(shapeErrorByEllipse.^2 + areaErrorByEllipse.^2));
        % We have to keep the fval non-infinite to keep BADS happy
        fval = min([fval realmax]);
    end


% Restore the warning state
warning(warningState);

% Assemble the sceneGeometry file to return
sceneGeometry = initialSceneGeometry;
sceneGeometry.cameraPosition.torsion = x(1);
sceneGeometry.cameraPosition.translation = x(2:4)';
sceneGeometry.eye.rotationCenters.azi = sceneGeometry.eye.rotationCenters.azi .* x(5) .* x(6);
sceneGeometry.eye.rotationCenters.ele = sceneGeometry.eye.rotationCenters.ele .* x(5) ./ x(6);
sceneGeometry.meta.estimateSceneParams.search.x = x';
sceneGeometry.meta.estimateSceneParams.search.options = options;
sceneGeometry.meta.estimateSceneParams.search.initialSceneGeometry = initialSceneGeometry;
sceneGeometry.meta.estimateSceneParams.search.ellipses = ellipses;
sceneGeometry.meta.estimateSceneParams.search.x0 = x0;
sceneGeometry.meta.estimateSceneParams.search.LB = LB;
sceneGeometry.meta.estimateSceneParams.search.UB = UB;
sceneGeometry.meta.estimateSceneParams.search.LBp = LBp;
sceneGeometry.meta.estimateSceneParams.search.UBp = UBp;
sceneGeometry.meta.estimateSceneParams.search.eyePoseLB = eyePoseLB;
sceneGeometry.meta.estimateSceneParams.search.eyePoseUB = eyePoseUB;
sceneGeometry.meta.estimateSceneParams.search.fVal = fVal;
sceneGeometry.meta.estimateSceneParams.search.centerDistanceErrorByEllipse = centerDistanceErrorByEllipse;
sceneGeometry.meta.estimateSceneParams.search.shapeErrorByEllipse = shapeErrorByEllipse;
sceneGeometry.meta.estimateSceneParams.search.areaErrorByEllipse = areaErrorByEllipse;
sceneGeometry.meta.estimateSceneParams.search.recoveredEyePoses = recoveredEyePoses;

end % local search function


function [] = saveSceneDiagnosticPlot(ellipses, Xedges, Yedges, eyePoseLB, eyePoseUB, sceneGeometry, sceneDiagnosticPlotFileName)
% Saves a plot that illustrates the sceneGeometry search results
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

% Prepare the figure
figHandle=figure('visible','off');
set(gcf,'PaperOrientation','landscape');

set(figHandle, 'Units','inches')
height = 6;
width = 11;

% The last two parameters of 'Position' define the figure size
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
    'eyePoseLB',eyePoseLB,'eyePoseUB',eyePoseUB),...
    1:1:size(ellipses,1),'UniformOutput',false);
projectedEllipses=vertcat(projectedEllipses{:});

% plot the projected ellipse centers
scatter(projectedEllipses(:,1),projectedEllipses(:,2),'o','filled', ...
    'MarkerFaceAlpha',2/8,'MarkerFaceColor',[0 0 1]);

% connect the centers with lines
ellipseRMSE = sceneGeometry.meta.estimateSceneParams.search.ellipseRMSE;
for ii=1:size(ellipses,1)
    lineAlpha = ellipseRMSE(ii)/max(ellipseRMSE);
    lineWeight = 0.5 + (ellipseRMSE(ii)/max(ellipseRMSE));
    ph=plot([projectedEllipses(ii,1) ellipses(ii,1)], ...
        [projectedEllipses(ii,2) ellipses(ii,2)], ...
        '-','Color',[1 0 0],'LineWidth', lineWeight);
    ph.Color(4) = lineAlpha;
end

% plot the estimated center of rotation of the eye
rotationCenterEllipse = pupilProjection_fwd([0 0 0 2], sceneGeometry);
plot(rotationCenterEllipse(1),rotationCenterEllipse(2), '+g', 'MarkerSize', 5);

% Calculate the plot limits
if ~isempty(Xedges)
    minX = min([projectedEllipses(:,1);ellipses(:,1);rotationCenterEllipse(1)]);
    maxX = max([projectedEllipses(:,1);ellipses(:,1);rotationCenterEllipse(1)]);
    minY = min([projectedEllipses(:,2);ellipses(:,2);rotationCenterEllipse(2)]);
    maxY = max([projectedEllipses(:,2);ellipses(:,2);rotationCenterEllipse(2)]);
    xPlotBounds = [(minX - (maxX-minX)/10) (maxX + (maxX-minX)/10) ];
    yPlotBounds = [(minY - (maxY-minY)/10) (maxY + (maxY-minY)/10) ];
else
    minX = min([projectedEllipses(:,1);ellipses(:,1);rotationCenterEllipse(1)]);
    maxX = max([projectedEllipses(:,1);ellipses(:,1);rotationCenterEllipse(1)]);
    minY = min([projectedEllipses(:,2);ellipses(:,2);rotationCenterEllipse(2)]);
    maxY = max([projectedEllipses(:,2);ellipses(:,2);rotationCenterEllipse(2)]);
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
shapeErrorVec = sceneGeometry.meta.estimateSceneParams.search.shapeErrorByEllipse;
shapeErrorVec = shapeErrorVec./sceneGeometry.constraintTolerance;
colorMatrix = zeros(3,size(ellipses,1));
colorMatrix(1,:)=1;
colorMatrix(2,:)= shapeErrorVec;
scatter(ellipses(:,1),ellipses(:,2),[],colorMatrix','o','filled');

% plot the estimated center of rotation of the eye
rotationCenterEllipse = pupilProjection_fwd([0 0 0 2], sceneGeometry);
plot(rotationCenterEllipse(1),rotationCenterEllipse(2), '+g', 'MarkerSize', 5);

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

% Add text to report the camera position parameters
xFinal = sceneGeometry.meta.estimateSceneParams.search.x;
myString = sprintf('torsion [deg] %4.1f; translation vector [mm] = %4.1f, %4.1f, %4.1f; rotation center scaling [joint, differential] = %4.2f, %4.2f',xFinal(1),xFinal(2),xFinal(3),xFinal(4),xFinal(5),xFinal(6));
text(0.5,1.0,myString,'Units','normalized','HorizontalAlignment','center')

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
areaErrorVec = sceneGeometry.meta.estimateSceneParams.search.areaErrorByEllipse;
areaErrorVec = abs(areaErrorVec)./sceneGeometry.constraintTolerance;
areaErrorVec = min([areaErrorVec ones(size(ellipses,1),1)],[],2);
colorMatrix = zeros(3,size(ellipses,1));
colorMatrix(1,:)=1;
colorMatrix(2,:)= areaErrorVec;
scatter(ellipses(:,1),ellipses(:,2),[],colorMatrix','o','filled');

% plot the estimated center of rotation of the eye
rotationCenterEllipse = pupilProjection_fwd([0 0 0 2], sceneGeometry);
plot(rotationCenterEllipse(1),rotationCenterEllipse(2), '+g', 'MarkerSize', 5);

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


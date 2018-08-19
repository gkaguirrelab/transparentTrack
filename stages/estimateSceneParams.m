function sceneGeometry = estimateSceneParams(pupilFileName, sceneGeometryFileName, varargin)
% Estimate camera position and eye rotation given image plane ellipses
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
%   minimizing the error in prediction of the center of those ellipses.
%
%   The search is conducted over 6 parameters, corresponding to torsional
%   camera rotation about the Z (depth) axis, three parameters of camera
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
%
% Optional key/value pairs (display and I/O):
%  'verbose'              - Logical. Default false.
%  'pupilFileToVideoSuffixSwitch' - Cell array that provides the suffix
%                           of the pupilData file and the suffix of the
%                           corresponding fit video file. This way, the fit
%                           video corresponding to the passed pupilData
%                           file can be found and used to create the
%                           ellipse array montage plot.
%
% Optional key/value pairs (flow control)
%  'useParallel'          - If set to true, use the MATLAB parallel pool
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
%                           separate worker if the parpool is available. If
%                           a value of zero is passed, then a sceneGeometry
%                           file and diagnostic plots are created using the
%                           midpoint of the passed bounds.
%  'nDiagnosticPlots'     - Scalar. Plots the n best solutions. If set to
%                           zero, no plots are saved.
%  'rankScaling'          - The relative influence of matching ellipse
%                           center, shape, and area (respectively) in
%                           judging the best solution.
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
p.addParameter('pupilFileToVideoSuffixSwitch',{'_pupil.mat','_gray.avi'},@iscell);

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
p.addParameter('fixationTargetArray',[],@(x)(isempty(x) | isnumeric(x)));
p.addParameter('nBinsPerDimension',4,@isnumeric);
p.addParameter('badFrameErrorThreshold',2, @isnumeric);
p.addParameter('nBADSsearches',10,@isnumeric);
p.addParameter('nDiagnosticPlots',5,@isnumeric);
p.addParameter('rankScaling',[4 2 1],@isnumeric);

% parse
p.parse(pupilFileName, sceneGeometryFileName, varargin{:})


%% Announce we are starting
if p.Results.verbose
    tic
    fprintf(['Estimating scene parameters. Started ' char(datetime('now')) '\n']);
end


%% Create initial sceneGeometry structure
initialSceneGeometry = createSceneGeometry(varargin{:});


%% Define the fixationTargetArray
fixationTargetArray = p.Results.fixationTargetArray;


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
    
    % If there is a fixationTargetArray, make sure it is the same length
    if ~isempty(fixationTargetArray)
        if length(fixationTargetArray) ~= length(ellipseArrayList)
            error('Unequal fixationTargetArray and ellipseArrayList');
        end
    end
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

    % If there is a fixationTargetArray, make sure it is empty
    if ~isempty(fixationTargetArray)
    	warning('Cannot use fixationTargetArray unless ellipseArrayList is defined');
    end
    fixationTargetArray=[];
end


%% Perform the search
% Inform the user
if p.Results.verbose
    fprintf(['Searching over scene geometry parameters.\n']);
    fprintf('| 0                      50                   100%% |\n');
    fprintf('.\n');
end

% If zero searches have been requested, then perform a fully bounded search
% with the solution parameters set to be the midpoint of the plausible
% bounds.
if p.Results.nBADSsearches==0
    sceneParams = (p.Results.sceneParamsLBp+p.Results.sceneParamsUBp)./2;
        sceneGeometry = ...
            performSceneSearch(initialSceneGeometry, ...
            ellipses(ellipseArrayList,:), ...
            ellipseFitRMSE(ellipseArrayList), ...
            fixationTargetArray, ...
            sceneParams, ...
            sceneParams, ...
            sceneParams, ...
            sceneParams, ...
            p.Results.eyePoseLB, ...
            p.Results.eyePoseUB);
        
        % Store the search results in the meta field
        tmpHold=sceneGeometry.meta.estimateSceneParams.search;
        sceneGeometry.meta.estimateSceneParams = p.Results;
        sceneGeometry.meta.estimateSceneParams.search = tmpHold;
        sceneGeometry.meta.estimateSceneParams.search.ellipseArrayList = ellipseArrayList';
        sceneGeometry.meta.estimateSceneParams.search.ellipseFitRMSE = ellipseFitRMSE(ellipseArrayList);
        rankOrder = 1;
        
        searchResults{1}=sceneGeometry;

else
    % Loop over the requested number of BADS searches
    searchResults = {};
    parfor (ss = 1:p.Results.nBADSsearches,nWorkers)
        
        searchResults{ss} = ...
            performSceneSearch(initialSceneGeometry, ...
            ellipses(ellipseArrayList,:), ...
            ellipseFitRMSE(ellipseArrayList), ...
            fixationTargetArray, ...
            p.Results.sceneParamsLB, ...
            p.Results.sceneParamsUB, ...
            p.Results.sceneParamsLBp, ...
            p.Results.sceneParamsUBp, ...
            p.Results.eyePoseLB, ...
            p.Results.eyePoseUB);
        
        % Save the ellipse details in the search results
        searchResults{ss}.meta.estimateSceneParams.search.ellipseArrayList = ellipseArrayList';
        searchResults{ss}.meta.estimateSceneParams.search.ellipseFitRMSE = ellipseFitRMSE(ellipseArrayList);

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
    
    % Obtain the "best" search result, which is defined as the smallest
    % product of the rankings of the three types of error (downweighting
    % the area error to have half the effect of shape and center)
    medianCenterErrorBySearch = cellfun(@(x) sqrt(nansum((x.meta.estimateSceneParams.search.centerDistanceErrorByEllipse).^2)),searchResults);
    medianShapeErrorBySearch = cellfun(@(x) sqrt(nansum((x.meta.estimateSceneParams.search.shapeErrorByEllipse).^2)),searchResults);
    medianAreaErrorBySearch = cellfun(@(x) sqrt(nansum((x.meta.estimateSceneParams.search.areaErrorByEllipse).^2)),searchResults);
    [~,centerErrorRank]  = ismember(medianCenterErrorBySearch,unique(medianCenterErrorBySearch));
    [~,shapeErrorRank]  = ismember(medianShapeErrorBySearch,unique(medianShapeErrorBySearch));
    [~,areaErrorRank]  = ismember(medianAreaErrorBySearch,unique(medianAreaErrorBySearch));
    b = p.Results.rankScaling;
    rankProduct = (centerErrorRank.*b(1)) .* ...
        (shapeErrorRank.*b(2)) .* ...
        (areaErrorRank.*b(3));
    [~,rankOrder]=sort(rankProduct);
    [~,bestSearchIdx] = min(rankProduct);
    sceneGeometry = searchResults{bestSearchIdx};

    % Store the search results in the meta field
    sceneGeometry.meta.estimateSceneParams = p.Results;
    sceneGeometry.meta.estimateSceneParams.allSearches = searchResults;
    sceneGeometry.meta.estimateSceneParams.bestSearchIdx = bestSearchIdx;
    sceneGeometry.meta.estimateSceneParams.search = searchResults{bestSearchIdx}.meta.estimateSceneParams.search;

end % Check for zero requested searches


%% Save the sceneGeometry file
if ~isempty(sceneGeometryFileName)
    save(sceneGeometryFileName,'sceneGeometry');
end

%% Save sceneGeometry diagnostics
if ~isempty(sceneGeometryFileName) && p.Results.nDiagnosticPlots~=0
    if p.Results.verbose
        fprintf('Creating diagnostic plots.\n');
    end
    [sceneGeomPath,sceneGeomName,] = fileparts(sceneGeometryFileName);
    diagnosticDirName = fullfile(sceneGeomPath,[sceneGeomName '_diagnostics']);
    if ~exist(diagnosticDirName, 'dir')
        mkdir(diagnosticDirName);
    else
        rmdir(diagnosticDirName, 's');
        mkdir(diagnosticDirName);
    end
    
    % Save the ellipse fit montage
    montageFileName = fullfile(diagnosticDirName,[sceneGeomName '_sceneDiagnosticMontage_ellipses.png']);
	grayVideoName = strrep(pupilFileName,p.Results.pupilFileToVideoSuffixSwitch{1},p.Results.pupilFileToVideoSuffixSwitch{2});    
    saveEllipseArrayMontage(sceneGeometry, ...
        ellipseArrayList, ...
        ellipses, ...
        grayVideoName, ...
        montageFileName ...
        );
        
    % Create a set of plots for the n best solutions
    parfor (ii = 1:min([length(rankOrder),p.Results.nDiagnosticPlots]))
        % Assemble the candidate sceneGeometry
        sceneDiagnosticPlotFileName = fullfile(diagnosticDirName,[sceneGeomName '_Rank' num2str(ii) '_Search' num2str(rankOrder(ii)) '_sceneDiagnosticPlot.pdf']);
        tmpSceneGeometry = searchResults{rankOrder(ii)};        
        
        % Create a sceneGeometry plot
        saveSceneDiagnosticPlot(...
            Xedges, Yedges,...
            p.Results.eyePoseLB, ...
            p.Results.eyePoseUB, ...
            tmpSceneGeometry,...
            sceneDiagnosticPlotFileName)
        
        % Create an eye model montage
        montageFileName = fullfile(diagnosticDirName,[sceneGeomName '_Rank' num2str(ii) '_Search' num2str(rankOrder(ii)) '_sceneDiagnosticMontage_eyeModel.png']);
        saveEyeModelMontage(tmpSceneGeometry, ...
            ellipseArrayList, ...
            ellipses, ...
            grayVideoName, ...
            montageFileName ...
            );
        
    end
end


%% alert the user that we are done with the routine
if p.Results.verbose
    toc
    fprintf('\n');
end


end % main function



%% LOCAL FUNCTIONS
function sceneGeometry = performSceneSearch(initialSceneGeometry, ellipses, ellipseFitRMSE, fixationTargetArray, LB, UB, LBp, UBp, eyePoseLB, eyePoseUB)
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

% Force the plausible bounds to be within the full bounds
UBp=min([UB';UBp'])';
LBp=max([LB';LBp'])';

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
warning('off','pupilProjection_fwd:ellipseFitFailed');

% Define nested variables for within the search
centerDistanceErrorByEllipse=zeros(size(ellipses,1),1);
shapeErrorByEllipse=zeros(size(ellipses,1),1);
areaErrorByEllipse=zeros(size(ellipses,1),1);
recoveredEyePoses =zeros(size(ellipses,1),4);
objEvalCounter = 0;
regParams = [];

% Tic
tic

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
        % Iterate the counter
        objEvalCounter = objEvalCounter+1;
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
                    'nMaxSearches',1);
        end
        % Objective function behavior varies depending upon if a fixation
        % target array list was provided
        if ~isempty(fixationTargetArray)
            % Compute the distance between the recovered eye rotation
            % angles and the fixation target positions
            [regParams,~,ErrorStats]=absor(...
                recoveredEyePoses(:,1:2)',...
                fixationTargetArray,...
                'weights',ellipseFitRMSE,...
                'doScale',false,...
                'doTrans',true);
            fval = ErrorStats.errlsq;
        else
            % Compute objective function as the RMSE of the distance
            % between the taget and modeled ellipses in shape and area
            fval = sqrt(mean(shapeErrorByEllipse.^2 + areaErrorByEllipse.^2));
            % We have to keep the fval non-infinite to keep BADS happy
            fval = min([fval realmax]);
        end
    end


% Restore the warning state
warning(warningState);

% Toc
searchTimeSecs = toc;

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
sceneGeometry.meta.estimateSceneParams.search.fixationTransform = regParams;
sceneGeometry.meta.estimateSceneParams.search.objEvalCounter = objEvalCounter;
sceneGeometry.meta.estimateSceneParams.search.searchTimeSecs = searchTimeSecs;

end % local search function


function [] = saveSceneDiagnosticPlot(Xedges, Yedges, eyePoseLB, eyePoseUB, sceneGeometry, sceneDiagnosticPlotFileName)
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

% Silence some errors that can arise during the forward projection
warningState = warning;
warning('off','pupilProjection_fwd:ellipseFitFailed');

% Obtain the set of ellipse parameters from the sceneGeometry structure
ellipses = sceneGeometry.meta.estimateSceneParams.search.ellipses;

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
ellipseFitRMSE = sceneGeometry.meta.estimateSceneParams.search.ellipseFitRMSE;
for ii=1:size(ellipses,1)
    lineAlpha = ellipseFitRMSE(ii)/max(ellipseFitRMSE);
    lineWeight = 0.5 + (ellipseFitRMSE(ii)/max(ellipseFitRMSE));
    if ~any(isnan(ellipses(ii,:))) && ~any(isnan(projectedEllipses(ii,:)))
        ph=plot([projectedEllipses(ii,1) ellipses(ii,1)], ...
            [projectedEllipses(ii,2) ellipses(ii,2)], ...
            '-','Color',[1 0 0],'LineWidth', lineWeight);
        ph.Color(4) = lineAlpha;
    end
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
else
    hold on
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
myString = sprintf('torsion [deg] = %4.1f; translation vector [mm] = %4.1f, %4.1f, %4.1f; rotation center scaling [joint, differential] = %4.2f, %4.2f',xFinal(1),xFinal(2),xFinal(3),xFinal(4),xFinal(5),xFinal(6));
text(0.5,1.0,myString,'Units','normalized','HorizontalAlignment','center')

% Add text to report the ellipse frames used
ellipseFrameList = num2str(sort(sceneGeometry.meta.estimateSceneParams.search.ellipseArrayList)');
myString = sprintf(['Ellipse frames (index from 1): [' ellipseFrameList ']']);
text(0.5,0.25,myString,'Units','normalized','HorizontalAlignment','center')


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
else
    hold on
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

% Restore the warning state
warning(warningState);

end % saveSceneDiagnosticPlot


function [] = saveEllipseArrayMontage(sceneGeometry, ellipseArrayList, allEllipses, grayVideoName, montageFileName)
% Saves a montage of the video frames illustrating the ellipses used for
% the sceneGeometry estimation.

% Sort the ellipse array list so that the frames appear in temporal order
ellipseArrayList = sort(ellipseArrayList);

% Check that the file exists
if exist(grayVideoName,'file') && ~isempty(ellipseArrayList)

    % Open the video object
    videoInObj = videoIOWrapper(grayVideoName,'ioAction','read');
    
    % Get the video properties
    videoSizeX = videoInObj.Width;
    videoSizeY = videoInObj.Height;
    nFrames = floor(videoInObj.Duration*videoInObj.FrameRate);
    
    % Define a variable to hold the selected frames
    framesToMontage = zeros(videoSizeY,videoSizeX,3,length(ellipseArrayList),'uint8');

    % Define a figure
    hFig = figure( 'Visible', 'off');
    hAxes = gca();
    
    % Loop through the frames 
    for ii = 1:length(ellipseArrayList)
        idx = ellipseArrayList(ii);
        videoInObj.CurrentTime = (idx - 1)/(videoInObj.FrameRate);
        sourceFrame = readFrame(videoInObj);
        sourceFrame = rgb2gray (sourceFrame);
        imshow(sourceFrame,'Border', 'tight','Parent',hAxes);
        hold on
        % Add the ellipse fit
        pFitImplicit = ellipse_ex2im(ellipse_transparent2ex(allEllipses(idx,:)));
        fh=@(x,y) pFitImplicit(1).*x.^2 +pFitImplicit(2).*x.*y +pFitImplicit(3).*y.^2 +pFitImplicit(4).*x +pFitImplicit(5).*y +pFitImplicit(6);
        fimplicit(fh,[1, videoSizeX, 1, videoSizeY],'Color', 'g','LineWidth',1);
        set(gca,'position',[0 0 1 1],'units','normalized')
        axis off;
        % Add blue grid lines
        lh = plot([1 videoSizeX],[round(videoSizeY/2) round(videoSizeY/2)],'-b','LineWidth',2);
        lh.Color=[0,0,1,0.5];
        lh = plot([round(videoSizeX/2) round(videoSizeX/2)],[1 videoSizeY],'-b','LineWidth',2);
        lh.Color=[0,0,1,0.5];
        % Get the frame
        drawnow;
        thisFrame=getframe(hFig);
        % Add a text label for the frame number
        frameLabel = sprintf('frame: %d',idx);
        thisFrame.cdata = insertText(thisFrame.cdata,[20 20],frameLabel,'FontSize',30);
        % Store the frame
        framesToMontage(:,:,:,ii) = thisFrame.cdata;
        % hold off
    end
    
    % Close the temporary figure
    close(hFig);
    
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
        'Color','w');

    % Turn off a warning that can occur during the montage step
    warningState = warning;
    warning('off','images:initSize:adjustingMag');
    
    % Create the montage
    montage(framesToMontage);
    
    % Restore the warning state
    warning(warningState);

    % Save the montage
    saveas(figHandle,montageFileName)
        
    % Close the figure
    close(figHandle)

    % Rotate the figure by 90 degrees clockwise, because I can't get the
    % MATLAB plotting routines to output the image how I want it.
    A = imread(montageFileName);
    A = rot90(A,3);
    imwrite(A,montageFileName);
    
    % close the video object
    clear videoInObj

    
end % There is a file to plot

end % saveEllipseArrayMontage



function [] = saveEyeModelMontage(sceneGeometry, ellipseArrayList, allEllipses, grayVideoName, montageFileName)
% Saves a montage with the model eye superimposed.

% Silence some errors that can arise during the forward projection
warningState = warning;
warning('off','pupilProjection_fwd:ellipseFitFailed');

% Sort the ellipse array list so that the frames appear in temporal order
ellipseArrayList = sort(ellipseArrayList);

% Check that the file exists
if exist(grayVideoName,'file') && ~isempty(ellipseArrayList)

    % Open the video object
    videoInObj = videoIOWrapper(grayVideoName,'ioAction','read');
    
    % Get the video properties
    videoSizeX = videoInObj.Width;
    videoSizeY = videoInObj.Height;
    nFrames = floor(videoInObj.Duration*videoInObj.FrameRate);
    
    % Define a variable to hold the selected frames
    framesToMontage = zeros(videoSizeY,videoSizeX,3,length(ellipseArrayList),'uint8');

    % Define a figure
    hFig = figure( 'Visible', 'off');
    hAxes = gca();
    
    % Loop through the frames and keep the matching ones
    for ii = 1:length(ellipseArrayList)
        idx = ellipseArrayList(ii);
        videoInObj.CurrentTime = (idx - 1)/(videoInObj.FrameRate);
        sourceFrame = readFrame(videoInObj);
        imshow(sourceFrame,'Border', 'tight','Parent',hAxes);
        hold on
        axis off;
        % Add the rendered eye model
        eyePose = sceneGeometry.meta.estimateSceneParams.search.recoveredEyePoses(ii,:);
        if ~any(isnan(eyePose))
            renderEyePose(eyePose, sceneGeometry, 'newFigure', false, ...
                'modelEyeLabelNames', {'retina' 'irisPerimeter' 'pupilEllipse' 'cornea'}, ...
                'modelEyePlotColors', {'.w' '.b' '-g' '.y'}, ...
                'modelEyeSymbolSizeScaler',1.5,...
                'modelEyeAlpha', 0.25);
        end
        % Get the frame
        drawnow;
        thisFrame=getframe(hFig);
        % Add a text label for the frame number
        frameLabel = sprintf('frame: %d',idx);
        thisFrame.cdata = insertText(thisFrame.cdata,[20 20],frameLabel,'FontSize',30);
        % Store the frame
        framesToMontage(:,:,:,ii) = thisFrame.cdata;
        % hold off
        hold off
    end

    % Close the temporary figure
    close(hFig);
    
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
        'Color','w');

    % Turn off a warning that can occur during the montage step
    warningState = warning;
    warning('off','images:initSize:adjustingMag');
    
    % Create the montage
    montage(framesToMontage);
    
    % Restore the warning state
    warning(warningState);

    % Save the montage
    saveas(figHandle,montageFileName)
        
    % Close the figure
    close(figHandle)

    % Rotate the figure by 90 degrees clockwise, because I can't get the
    % MATLAB plotting routines to output the image how I want it.
    A = imread(montageFileName);
    A = rot90(A,3);
    imwrite(A,montageFileName);
    
    % close the video object
    clear videoInObj
    
end % There is a file to plot

% Restore the warning state
warning(warningState);

end % saveEyeModelMontage

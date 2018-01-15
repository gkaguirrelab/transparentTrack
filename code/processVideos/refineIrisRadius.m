function sceneGeometry = refineIrisRadius(pupilFileName, sceneGeometryFileName, varargin)
% Adjust the model iris radius to best fit the image plane data 
%
% Description:
%   This function searches over values for the radius of the iris in the
%   scene to best fit the border of the iris on the image plane. The
%   routine takes as input a pupilData and sceneGeometry file. A set of
%   frames from the pupilData file are selected (or passed). The iris that
%   best identifies an image intensity boundary between the iris and sclera
%   is found. This value is placed into the sceneGeometry structure and
%   saved.
%
% Inputs:
%	pupilFileName         - Full path to a pupilData file, or a cell array
%                           of such paths. If a cell array is provided, the
%                           ellipse data from each pupilData file will be
%                           loaded and concatenated. If set to empty, a
%                           sceneGeometry structure with default values
%                           will be returned.
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
%  'irisRadiusLB'         - Scalar
%  'irisRadiusUB'         - Scalar
%  'whichEllipseFitField' - Identifies the field in pupilData that contains
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
p.addRequired('pupilFileName',@(x)(isempty(x) | iscell(x) | ischar(x)));
p.addRequired('sceneGeometryFileName',@(x)(isempty(x) | ischar(x)));

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
p.addParameter('irisRadiusLB',12.5,@isnumeric);
p.addParameter('irisRadiusUB',15.0,@isnumeric);
p.addParameter('ellipseFitLabel', 'radiusSmoothed', @ischar);
p.addParameter('ellipseArrayList',[],@(x)(isempty(x) | isnumeric(x)));
p.addParameter('nBinsPerDimension',10,@isnumeric);

% parse
p.parse(pupilFileName, sceneGeometryFileName, varargin{:})



%% Announce we are starting
if strcmp(p.Results.verbosity,'full')
    tic
    fprintf(['Refining iris readius. Started ' char(datetime('now')) '\n']);
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


%% Load pupil data
if iscell(pupilFileName)
    ellipses = [];
    ellipseFitSEM = [];
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


%% Identify the ellipses that will guide the sceneGeometry estimation
% If not supplied, we will generate a list of ellipses to use for the
% estimation.
if ~isempty(p.Results.ellipseArrayList)
    ellipseArrayList = p.Results.ellipseArrayList;
    Xedges = [];
    Yedges = [];
else
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

%% Load the sceneGeometry
dataLoad = load(p.Results.sceneGeometryFileName);
initialSceneGeometry = dataLoad.sceneGeometry;
clear dataLoad


%% Perform the search
% Call out to the local function that performs the serach
sceneGeometry = ...
    performIrisRadiusSearch(initialSceneGeometry, ...
    ellipses(ellipseArrayList,:), ...
    errorWeights, ...
    p.Results.irisRadiusLB, ...
    p.Results.irisRadiusUB);

% add additional search and meta field info to sceneGeometry
sceneGeometry.meta.refineIris.parameters = p.Results;
sceneGeometry.meta.refineIris.search.ellipseArrayList = ellipseArrayList';


%% Save the sceneGeometry file
if ~isempty(sceneGeometryFileName)
    save(sceneGeometryFileName,'sceneGeometry');
end


%% alert the user that we are done with the routine
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

function sceneGeometry = performIrisRadiusSearch(initialSceneGeometry, ellipses, errorWeights, irisRadiusLB, irisRadiusUB)
% Search for best fitting irisRadius
%
% Description:
%
%

% Set the error form
errorForm = 'SSE';

% Extract the initial search point from initialSceneGeometry
x0 = initialSceneGeometry.eye.irisRadius;

% Define search options
options = optimoptions(@patternsearch, ...
    'Display','off',...
    'AccelerateMesh',false,...
    'Cache','on',...
    'CompleteSearch','on',...
    'UseParallel', true, ...
    'FunctionTolerance',1e-6);

% Define anonymous functions for the objective and constraint
objectiveFun = @objfun; % the objective function, nested below

% Define nested variables for within the search
centerDistanceErrorByEllipse=[];
shapeErrorByEllipse=[];
areaErrorByEllipse=[];

[x, fVal] = patternsearch(objectiveFun, x0,[],[],[],[],sceneParamsLB,sceneParamsUB,[],options);
    function fval = objfun(x)
        candidateSceneGeometry = initialSceneGeometry;
        candidateSceneGeometry.extrinsicTranslationVector = x(1:3);
        candidateSceneGeometry.eye.centerOfRotation(1) = -x(4);
        [~, ~, centerDistanceErrorByEllipse, shapeErrorByEllipse, areaErrorByEllipse] = ...
            arrayfun(@(x) pupilProjection_inv...
            (...
                ellipses(x,:),...
                candidateSceneGeometry,...
                'constraintTolerance', candidateSceneGeometry.constraintTolerance,...
                'eyeParamsLB',eyeParamsLB,...
                'eyeParamsUB',eyeParamsUB...
            ),...
            1:1:size(ellipses,1),'UniformOutput',false);
        
        % Now compute objective function as the RMSE of the distance
        % between the taget and modeled ellipses
        centerDistanceErrorByEllipse = cell2mat(centerDistanceErrorByEllipse)';
        shapeErrorByEllipse = cell2mat(shapeErrorByEllipse)';
        areaErrorByEllipse = cell2mat(areaErrorByEllipse)';
        switch errorForm
            case 'SSE'
                fval=sum((centerDistanceErrorByEllipse.*(shapeErrorByEllipse.*100+1).*(areaErrorByEllipse.*100+1).*errorWeights).^2);
            case 'RMSE'
                fval = mean((centerDistanceErrorByEllipse.*(shapeErrorByEllipse.*100+1).*(areaErrorByEllipse.*100+1).*errorWeights).^2)^(1/2);
            otherwise
                error('I do not recognize that error form');
        end
        
    end

% Assemble the sceneGeometry file to return
sceneGeometry.radialDistortionVector = initialSceneGeometry.radialDistortionVector;
sceneGeometry.intrinsicCameraMatrix = initialSceneGeometry.intrinsicCameraMatrix;
sceneGeometry.extrinsicTranslationVector = x(1:3);
sceneGeometry.extrinsicRotationMatrix = initialSceneGeometry.extrinsicRotationMatrix;
sceneGeometry.constraintTolerance = initialSceneGeometry.constraintTolerance;
sceneGeometry.eye = initialSceneGeometry.eye;
sceneGeometry.eye.centerOfRotation(1) = -x(4);
sceneGeometry.meta.refineIris.search.options = options;
sceneGeometry.meta.refineIris.search.errorForm = errorForm;
sceneGeometry.meta.refineIris.search.initialSceneGeometry = initialSceneGeometry;
sceneGeometry.meta.refineIris.search.ellipses = ellipses;
sceneGeometry.meta.refineIris.search.errorWeights = errorWeights;
sceneGeometry.meta.refineIris.search.sceneParamsLB = sceneParamsLB;
sceneGeometry.meta.refineIris.search.sceneParamsUB = sceneParamsUB;
sceneGeometry.meta.refineIris.search.eyeParamsLB = eyeParamsLB;
sceneGeometry.meta.refineIris.search.eyeParamsUB = eyeParamsUB;
sceneGeometry.meta.refineIris.search.fVal = fVal;
sceneGeometry.meta.refineIris.search.centerDistanceErrorByEllipse = centerDistanceErrorByEllipse;
sceneGeometry.meta.refineIris.search.shapeErrorByEllipse = shapeErrorByEllipse;
sceneGeometry.meta.refineIris.search.areaErrorByEllipse = areaErrorByEllipse;

end % local search function



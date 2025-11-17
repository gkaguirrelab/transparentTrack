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
%   coordinates. Additional information regarding the MATLAB image
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
%  'verbose'              - Logical. Default false.
%
% Optional key/value pairs (flow control)
%  'nFrames'              - Analyze fewer than the total number of frames.
%  'startFrame'           - Which frame to start on
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
%  'ellipseTransparentLB/UB' - Define the hard upper and lower boundaries
%                           for the ellipse fit, in units of pixels of the
%                           video. The default values selected here
%                           represent the physical and mathematical limits,
%                           as the constraint for the fit will be provided
%                           by the scene geometry. A mild constraint (0.6)
%                           is placed upon the eccentricity, corresponding
%                           to an aspect ration of 4:5.
%  'eyePoseLB'            - Lower bound on the eyePose
%  'eyePoseUB'            - Upper bound on the eyePose
%  'glintFileName'        - Full path to a glint file. When available, the
%                           glint is used to constrain the eyePose that is
%                           found to fit the pupil perimeter.
%  'sceneGeometryFileName' - Full path to a sceneGeometry file. When the
%                           sceneGeometry is available, fitting is
%                           performed in terms of eye parameters instead of
%                           ellipse parameters
%  'fitLabel'             - The field name in the pupilData structure where
%                           the results of the fitting will be stored.
%  'relativeCameraPositionFileName' - Char. This is the full path to a
%                           relativeCameraPosition.mat file that provides
%                           the position of the camera at each video frame
%                           relative to the initial position of the camera.
%
% Outputs:
%	pupilData             - A structure with multiple fields corresponding
%                           to the parameters and errors of the fit.
%                           Different field names are used depending upon
%                           if a sceneGeometry constraint was or was not
%                           used.
%

%% Parse vargin for options passed here
p = inputParser; p.KeepUnmatched = true; p.PartialMatching = false;

% Required
p.addRequired('perimeterFileName',@ischar);
p.addRequired('pupilFileName',@ischar);

% Optional display and I/O params
p.addParameter('verbose',false,@islogical);

% Optional flow control params
p.addParameter('nFrames',Inf,@isnumeric);
p.addParameter('startFrame',1,@isnumeric);
p.addParameter('useParallel',false,@islogical);
p.addParameter('nWorkers',[],@(x)(isempty(x) | isnumeric(x)));

% Optional environment params
p.addParameter('tbSnapshot',[],@(x)(isempty(x) | isstruct(x)));
p.addParameter('timestamp',char(datetime('now')),@ischar);
p.addParameter('hostname',char(java.lang.System.getProperty('user.name')),@ischar);

% Optional analysis params
p.addParameter('ellipseTransparentLB',[0,0,800,0,0],@(x)(isempty(x) || isnumeric(x)));
p.addParameter('ellipseTransparentUB',[640,480,20000,0.6,pi],@(x)(isempty(x) || isnumeric(x)));
p.addParameter('eyePoseLB',[-89,-89,0,0.1],@isnumeric);
p.addParameter('eyePoseUB',[89,89,0,5],@isnumeric);
p.addParameter('cameraTransBounds',[5; 5; 0],@isnumeric);
p.addParameter('sceneGeometryFileName',[],@(x)(isempty(x) || ischar(x)));
p.addParameter('confidenceThreshold',0.75,@isnumeric);
p.addParameter('glintFileName',[],@(x)(isempty(x) || ischar(x)));
p.addParameter('glintTol',1,@isnumeric);
p.addParameter('fitLabel',[],@(x)(isempty(x) | ischar(x)));
p.addParameter('relativeCameraPositionFileName',[],@ischar);


%% Parse and check the parameters
p.parse(perimeterFileName, pupilFileName, varargin{:});

nEllipseParams=5; % 5 params in the transparent ellipse form
nEyePoseParams=4; % [azimuth, elevation, torsion, radius]
nHeadTransParams=3; % [horizontal; vertical; depth]


%% Load data
% Load the pupil perimeter data. It will be a structure variable
% "perimeter", with the fields .data and .meta
load(perimeterFileName,'perimeter');

% Load the pupilData file if it exists and we have been given a
% sceneGeometry file. In this circumstance, we are loading the file that
% contains the initial fit, and we are not going to add the scene
% constrained fit.
if exist(p.Results.pupilFileName, 'file')==2 && ~isempty(p.Results.sceneGeometryFileName)
    load(pupilFileName,'pupilData');
else
    pupilData=[];
end

% Load the sceneGeometry file if passed
if ~isempty(p.Results.sceneGeometryFileName)
    % Load the sceneGeometry file
    load(p.Results.sceneGeometryFileName,'sceneGeometry');
else
    sceneGeometry = [];
end

% Load the glint file if passed
glintTol = p.Results.glintTol;
if ~isempty(p.Results.glintFileName)
    % Load the glintData file
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

% Load the relativeCameraPosition file if passed and it exists
if ~isempty(p.Results.relativeCameraPositionFileName)
    if exist(p.Results.relativeCameraPositionFileName, 'file')==2
        load(p.Results.relativeCameraPositionFileName,'relativeCameraPosition');
    else
        relativeCameraPosition=[];
    end
else
    relativeCameraPosition=[];
end

% Select the appropriate field of relativeCameraTransition, or synthesize
% one if not available
if ~isempty(relativeCameraPosition)
    cameraTransVec = ...
        relativeCameraPosition.(relativeCameraPosition.currentField).values;
else
    cameraTransVec = zeros(nHeadTransParams,size(perimeter.data,1)+1);
end

% determine how many frames we will process
startFrame = p.Results.startFrame;
nFrames = size(perimeter.data,1);


%% Establish a label for this analysis stage
if isempty(p.Results.fitLabel)
    if isempty(sceneGeometry)
        fitLabel = 'initial';
    else
        fitLabel = 'sceneConstrained';
    end
else
    fitLabel = p.Results.fitLabel;
end


%% Set up the parallel pool
if p.Results.useParallel
    nWorkers = startParpool( p.Results.nWorkers, p.Results.verbose );
else
    nWorkers=0;
end


%% Prepare for the loop

% Recast perimeter into a sliced cell array to reduce parfor broadcast
% overhead
frameCellArray = perimeter.data;

% Set-up other variables to be non-broadcast
verbose = p.Results.verbose;
ellipseTransparentLB = p.Results.ellipseTransparentLB;
ellipseTransparentUB = p.Results.ellipseTransparentUB;
eyePoseLB = p.Results.eyePoseLB;
eyePoseUB = p.Results.eyePoseUB;
cameraTransBounds = p.Results.cameraTransBounds;
confidenceThreshold = p.Results.confidenceThreshold;

% Alert the user
if p.Results.verbose
    tic
    fprintf(['Ellipse fitting to pupil perimeter. Started ' char(datetime('now')) '\n']);
    fprintf('| 0                      50                   100%% |\n');
    fprintf('.\n');
end

% Store the warning state
warnState = warning();

%% Loop through the frames
parfor (ii = startFrame:nFrames, nWorkers)
%for ii = startFrame:nFrames
%for ii = 16341:16345

    % Update progress
    if verbose
        if mod(ii,round((nFrames-startFrame)/50))==0
            fprintf('\b.\n');
        end
    end
    
    % Initialize the results variables
    ellipseParamsTransparent=NaN(1,nEllipseParams);
    objectiveError=NaN(1);
    cameraTrans=NaN(nHeadTransParams,1);
    eyePose=NaN(1,nEyePoseParams);
    fitAtBound=false;

    % Get the camera translation for this frame
    cameraTransInitial = cameraTransVec(:,ii);

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

    % Ensure that the perimeter points are a column vector
    if ~iscolumn(Xp)
        Xp = Xp'; Yp = Yp';
    end

    % fit an ellipse to the boundary (if any points exist)
    if length(Xp) > 1 && length(Yp) > 1
        
        % Turn off expected warnings
        warning('off','projectModelEye:ellipseFitFailed');
        warning('off','gkaModelEye:pupilEllipseFit');
        warning('off','MATLAB:nearlySingularMatrix');
        warning('off','MATLAB:singularMatrix');
        
        % Fit approach depends upon whether or not we have sceneGeometry
        if isempty(sceneGeometry)
            % No sceneGeometry. Fit an ellipse to the perimeter.
            [ellipseParamsTransparent, objectiveError, ~, fitAtBound] = ...
                constrainedEllipseFit(Xp, Yp, ...
                ellipseTransparentLB, ...
                ellipseTransparentUB, ...
                []);
        else
            
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
            
            % Find the eyePose parameters that best fit the pupil
            % perimeter. This can take a few seconds.
            [eyePose, cameraTrans, objectiveError, ellipseParamsTransparent, fitAtBound] = ...
                eyePoseEllipseFit(Xp, Yp, glintCoord, sceneGeometry, ...
                'glintTol',glintTol, 'cameraTransX0',cameraTransInitial,...
                'cameraTransBounds',thisFrameCameraTransBounds,...
                'eyePoseLB', eyePoseLB, 'eyePoseUB', eyePoseUB);

        end
        
        % Restore the warning state
        warning(warnState);
        
    end % check if there are pupil boundary data to be fit
    
    % store results
    loopVar_ellipseParamsTransparent(ii,:) = ellipseParamsTransparent';
    loopVar_objectiveError(ii) = objectiveError;
    loopVar_fitAtBound(ii) = fitAtBound;
    if ~isempty(sceneGeometry)
        loopVar_eyePoses(ii,:) = eyePose';
        loopVar_cameraTransVec(:,ii) = cameraTrans;
    end
    
    
end % loop over frames

% alert the user that we are done with the fit loop
if p.Results.verbose
    toc
    fprintf('\n');
end

%% Clean up and save

% Clear out any old results in this fit label field
pupilData.(fitLabel) = [];

% Store the identity of the most recently produced field of data
pupilData.currentField = fitLabel;

% Store the ellipse fit data in informative fields
pupilData.(fitLabel).ellipses.values = loopVar_ellipseParamsTransparent;
pupilData.(fitLabel).ellipses.RMSE = loopVar_objectiveError';
if isempty(sceneGeometry)
    pupilData.(fitLabel).ellipses.fitAtBound = loopVar_fitAtBound';
else
    pupilData.(fitLabel).eyePoses.fitAtBound = loopVar_fitAtBound';
end
pupilData.(fitLabel).ellipses.meta.ellipseForm = 'transparent';
pupilData.(fitLabel).ellipses.meta.labels = {'x','y','area','eccentricity','theta'};
pupilData.(fitLabel).ellipses.meta.units = {'pixels','pixels','squared pixels','non-linear eccentricity','rads'};
pupilData.(fitLabel).ellipses.meta.coordinateSystem = 'intrinsic image';
if ~isempty(sceneGeometry)
    % Update the pupilData
    pupilData.(fitLabel).eyePoses.values = loopVar_eyePoses;
    pupilData.(fitLabel).eyePoses.meta.labels = {'azimuth','elevation','torsion','pupil radius'};
    pupilData.(fitLabel).eyePoses.meta.units = {'deg','deg','deg','mm'};
    pupilData.(fitLabel).eyePoses.meta.coordinateSystem = 'head fixed (extrinsic)';
    
    % Update the relativeCameraPosition
    relativeCameraPosition.(fitLabel).values = loopVar_cameraTransVec;
    relativeCameraPosition.(fitLabel).meta = p.Results;
    relativeCameraPosition.currentField = fitLabel;
end

% If the perimeter variable has an "instructions" field, copy it over to
% pupilData
if isfield(perimeter,'instructions')
    pupilData.instructions = perimeter.instructions;
end

% add meta data
pupilData.(fitLabel).meta = p.Results;

% save the ellipse fit results
save(p.Results.pupilFileName,'pupilData')

% save the relativeCameraPosition
if ~isempty(sceneGeometry) && ~isempty(p.Results.relativeCameraPositionFileName)
    save(p.Results.relativeCameraPositionFileName,'relativeCameraPosition')
end

end % function


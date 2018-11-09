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
%   coordinates. Additinal information regarding the MATLAB image
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
%  'sceneGeometryFileName' - Full path to a sceneGeometry file. When the
%                           sceneGeometry is available, fitting is
%                           performed in terms of eye parameters instead of
%                           ellipse parameters
%  'badFrameErrorThreshold' - This values is passed to eyePoseEllipseFit(),
%                           which is the function used when performing
%                           scene constrained fitting. The value is used to
%                           detect a possible local minimum when performing
%                           the eye pose search.
%  'fitLabel'             - The field name in the pupilData structure where
%                           the results of the fitting will be stored.
%
% Outputs:
%	pupilData             - A structure with multiple fields corresponding
%                           to the parameters and errors of the fit.
%                           Different field names are used depending upon
%                           if a sceneGeometry constraint was or was not
%                           used.
%

%% Parse vargin for options passed here
p = inputParser; p.KeepUnmatched = true;

% Required
p.addRequired('perimeterFileName',@ischar);
p.addRequired('pupilFileName',@ischar);

% Optional display and I/O params
p.addParameter('verbose',false,@islogical);

% Optional flow control params
p.addParameter('nFrames',Inf,@isnumeric);
p.addParameter('useParallel',false,@islogical);
p.addParameter('nWorkers',[],@(x)(isempty(x) | isnumeric(x)));

% Optional environment params
p.addParameter('tbSnapshot',[],@(x)(isempty(x) | isstruct(x)));
p.addParameter('timestamp',char(datetime('now')),@ischar);
p.addParameter('hostname',char(java.lang.System.getProperty('user.name')),@ischar);
p.addParameter('username',char(java.net.InetAddress.getLocalHost.getHostName),@ischar);

% Optional analysis params
p.addParameter('ellipseTransparentLB',[0,0,800,0,0],@(x)(isempty(x) || isnumeric(x)));
p.addParameter('ellipseTransparentUB',[640,480,20000,0.6,pi],@(x)(isempty(x) || isnumeric(x)));
p.addParameter('eyePoseLB',[-89,-89,0,0.1],@isnumeric);
p.addParameter('eyePoseUB',[89,89,0,5],@isnumeric);
p.addParameter('sceneGeometryFileName',[],@(x)(isempty(x) || ischar(x)));
p.addParameter('badFrameErrorThreshold',2, @isnumeric);
p.addParameter('fitLabel',[],@(x)(isempty(x) | ischar(x)));
p.addParameter('relativeCameraPositionFileName',[],@ischar);


%% Parse and check the parameters
p.parse(perimeterFileName, pupilFileName, varargin{:});

nEllipseParams=5; % 5 params in the transparent ellipse form
nEyePoseParams=4; % 4 eyePose values (azimuth, elevation, torsion, radius) 


%% Load data
% Load the pupil perimeter data. It will be a structure variable
% "perimeter", with the fields .data and .meta
dataLoad=load(perimeterFileName);
perimeter=dataLoad.perimeter;
clear dataLoad

% Load the pupilData file if it exists and we have been given a
% sceneGeometry file. In this circumstance, we are loading the file that
% contains the initial fit, and we are not going to add the scene
% constrained fit.
if exist(p.Results.pupilFileName, 'file')==2 && ~isempty(p.Results.sceneGeometryFileName)
    dataLoad=load(pupilFileName);
    pupilData=dataLoad.pupilData;
    clear dataLoad
else
    pupilData=[];
end

% Load the relativeCameraPosition file if passed and it exists
if ~isempty(p.Results.relativeCameraPositionFileName)
    if exist(p.Results.pupilFileName, 'file')==2
        dataLoad=load(p.Results.relativeCameraPositionFileName);
        relativeCameraPosition=dataLoad.relativeCameraPosition;
        clear dataLoad
    else
        relativeCameraPosition=[];
    end
else
    relativeCameraPosition=[];
end

% determine how many frames we will process
if p.Results.nFrames == Inf
    nFrames=size(perimeter.data,1);
else
    nFrames = p.Results.nFrames;
end


%% Set up the parallel pool
if p.Results.useParallel
    nWorkers = startParpool( p.Results.nWorkers, p.Results.verbose );
else
    nWorkers=0;
end

% Load the sceneGeometry file if passed
if ~isempty(p.Results.sceneGeometryFileName)
    % Load the sceneGeometry file
    dataLoad=load(p.Results.sceneGeometryFileName);
    sceneGeometry=dataLoad.sceneGeometry;
    clear dataLoad
    % An earlier version of the code defined a non-zero iris thickness. We
    % force this to zero here to speed computation
    sceneGeometry.eye.iris.thickness=0;
else
    sceneGeometry = [];
end


%% Calculate an ellipse fit for each video frame

% Recast perimeter into a sliced cell array to reduce parfor broadcast
% overhead
frameCellArray = perimeter.data(1:nFrames);

% Set-up other variables to be non-broadcast
verbose = p.Results.verbose;
ellipseTransparentLB = p.Results.ellipseTransparentLB;
ellipseTransparentUB = p.Results.ellipseTransparentUB;
eyePoseLB = p.Results.eyePoseLB;
eyePoseUB = p.Results.eyePoseUB;
badFrameErrorThreshold = p.Results.badFrameErrorThreshold;

% Alert the user
if p.Results.verbose
    tic
    fprintf(['Ellipse fitting to pupil perimeter. Started ' char(datetime('now')) '\n']);
    fprintf('| 0                      50                   100%% |\n');
    fprintf('.\n');
end

% Store the warning state
warnState = warning();

% Loop through the frames
parfor (ii = 1:nFrames, nWorkers)
    %for ii = 1:nFrames

    % Update progress
    if verbose
        if mod(ii,round(nFrames/50))==0
            fprintf('\b.\n');
        end
    end
    
    % Initialize the results variables
    ellipseParamsTransparent=NaN(1,nEllipseParams);
    objectiveError=NaN(1);
    eyePose=NaN(1,nEyePoseParams);
    fitAtBound=NaN(1);
    
    % get the boundary points
    Xp = frameCellArray{ii}.Xp;
    Yp = frameCellArray{ii}.Yp;
    
    % fit an ellipse to the boundary (if any points exist)
    if ~isempty(Xp) && ~isempty(Yp)

        % Turn off expected warnings
        warning('off','pupilProjection_fwd:ellipseFitFailed');
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
            % We do have sceneGeometry. 
            adjustedSceneGeometry = sceneGeometry;
            % If a relativeCameraPosition is defined, update the
            % sceneGeometry
            if ~isempty(relativeCameraPosition)
                cameraPosition = sceneGeometry.cameraPosition.translation;
                cameraPosition = cameraPosition + relativeCameraPosition(:,ii);
                adjustedSceneGeometry.cameraPosition.translation = cameraPosition;
            end
            % Find the eyePose parameters that best fit the pupil perimeter
            [eyePose, objectiveError, ellipseParamsTransparent, fitAtBound] = ...
                eyePoseEllipseFit(Xp, Yp, adjustedSceneGeometry, ...
                'eyePoseLB', eyePoseLB, 'eyePoseUB', eyePoseUB, ...
                'repeatSearchThresh', badFrameErrorThreshold);
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
    end
    
end % loop over frames

% alert the user that we are done with the fit loop
if p.Results.verbose
    toc
    fprintf('\n');
end

%% Clean up and save

% Establish a label to save the fields of the ellipse fit data
if isempty(p.Results.fitLabel)
    if isempty(sceneGeometry)
        fitLabel = 'initial';
    else
        fitLabel = 'sceneConstrained';
    end
else
    fitLabel = p.Results.fitLabel;
end

% Clear out any old results in this fit label field
pupilData.(fitLabel) = [];

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
    pupilData.(fitLabel).eyePoses.values = loopVar_eyePoses;
    pupilData.(fitLabel).eyePoses.meta.labels = {'azimuth','elevation','torsion','pupil radius'};
    pupilData.(fitLabel).eyePoses.meta.units = {'deg','deg','deg','mm'};
    pupilData.(fitLabel).eyePoses.meta.coordinateSystem = 'head fixed (extrinsic)';
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

end % function


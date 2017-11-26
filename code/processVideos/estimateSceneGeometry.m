function sceneGeometry = estimateSceneGeometry(pupilFileName, sceneGeometryFileName, varargin)
% sceneGeometry = estimateSceneGeometry(pupilFileName, sceneGeometryFileName)
%
% This function searches over the set of ellipses in the passed pupil file
% to estimate the sceneGeometry features in units of pixels on the scene.
% The routine identifies the [X, Y, Z] coordinates on the scene plane of a
% the center of rotation of an eye that would minimize the center of
% projection implied by each ellipse in the pupil file.
%
% Different projection models can be used to guide this calculation: The
% orthogonal model assumes that ellipses  on the scene plane are orthogonal
% projections of a circular pupil the center of which rotates around the
% eye center.
%
% The x, y coordinates on the scene plane of the center of the most
% circular (least eccentric) ellipse are taken as the initial guess for the
% center of rotation of the eye.
%
%
% Input
%	pupilFileName - Full path to a pupilData file
%   sceneGeometryFileName -  Full path to the file in which the
%      sceneGeometry data should be saved.
%
% Output
%	sceneGeometry - A structure with the fields
%       eyeCenter.X - X coordinate of the eye center (i.e. the assumed
%           center of rotation of the pupil) on the scene plane.
%       eyeCenter.Y - Y coordinate of the eye center (i.e. the assumed
%           center of rotation of the pupil) on the scene plane.
%       eyeCenter.Z - the orthogonal distance for the eye center from the
%           scene plane.
%       eyeRadius - radius of the eye in pixels
%       meta - information regarding the analysis, including units.
%
% Options (analysis)
%   projectionModel - options include 'orthogonal' and 'perspective'
%   eyeRadius - an assigned valye for eye radius. When operating
%       with the orthogonal projection model, this value is stored and
%       used for all subsequent calculations. When operating in the
%       perspective model, the eyeRadius is estimated from the data.
%
% Options (verbosity and display)
%   verbosity - controls console status updates
%   displayMode - when set to true, displays the results of the fit
%
% Optional key/value pairs (flow control)
%  'nFrames' - analyze fewer than the total number of frames.
%  'startFrame' - which frame to start on
%
% Options (environment)
%   tbSnapshot - the passed tbSnapshot output that is to be saved along
%      with the data
%   timestamp / username / hostname - these are automatically derived and
%      saved within the p.Results structure.
%


%% input parser
p = inputParser; p.KeepUnmatched = true;

% required input
p.addRequired('pupilFileName',@isstr);
p.addRequired('sceneGeometryFileName',@isstr);

% Optional analysis params
p.addParameter('projectionModel','pseudoPerspective',@ischar);
p.addParameter('eyeRadius',125,@isnumeric);
p.addParameter('cameraDistanceInPixels',1200,@isnumeric);
p.addParameter('sceneGeometryLowerBounds',[0, 0, 1325, 25],@isnumeric);
p.addParameter('sceneGeometryUpperBounds',[640, 480, 1325, 500],@isnumeric);
p.addParameter('numberFitDivisions',2,@isnumeric);
p.addParameter('whichFitFieldMean','ellipseParamsUnconstrained_mean',@ischar);
p.addParameter('whichFitFieldError','ellipseParamsUnconstrained_rmse',@ischar);

% verbosity and plotting control
p.addParameter('verbosity', 'none', @isstr);
p.addParameter('sceneDiagnosticPlotFileName', [],@(x)(isempty(x) | ischar(x)));
p.addParameter('sceneDiagnosticPlotSizeXY', [640 480],@isnumeric);

% flow control
p.addParameter('nFrames',Inf,@isnumeric);
p.addParameter('startFrame',1,@isnumeric);

% Environment parameters
p.addParameter('tbSnapshot',[],@(x)(isempty(x) | isstruct(x)));
p.addParameter('timestamp',char(datetime('now')),@ischar);
p.addParameter('username',char(java.lang.System.getProperty('user.name')),@ischar);
p.addParameter('hostname',char(java.net.InetAddress.getLocalHost.getHostName),@ischar);

% parse
p.parse(pupilFileName, sceneGeometryFileName, varargin{:})

%% main

%% Announce we are starting
if strcmp(p.Results.verbosity,'full')
    tic
    fprintf(['Estimating scene geometry from pupil ellipses. Started ' char(datetime('now')) '\n']);
end

% load pupil data
load(pupilFileName)
if p.Results.nFrames ~= Inf
    ellipses = pupilData.(p.Results.whichFitFieldMean)(p.Results.startFrame:p.Results.nFrames,:);
else
    ellipses = pupilData.(p.Results.whichFitFieldMean)(p.Results.startFrame:end,:);
end

% find the most circular ellipse and use the X Y coordinate as the initial
% guess for the CoP XY coordinates
[~, minEccentricityIdx] = min(pupilData.(p.Results.whichFitFieldMean)(:,4));

% Identify the [X Y] coordinates of the most circular ellipse. Z is set to
% an initial value of zero
switch p.Results.projectionModel
    case 'orthogonal'
        x0 = [pupilData.(p.Results.whichFitFieldMean)(minEccentricityIdx,1) ...
            pupilData.(p.Results.whichFitFieldMean)(minEccentricityIdx,2) ...
            0 ...
            p.Results.eyeRadius];
    case 'pseudoPerspective'
        x0 = [pupilData.(p.Results.whichFitFieldMean)(minEccentricityIdx,1) ...
            pupilData.(p.Results.whichFitFieldMean)(minEccentricityIdx,2) ...
            p.Results.cameraDistanceInPixels + p.Results.eyeRadius ...
            p.Results.eyeRadius];
    otherwise
        error('I do not recognize that perspective model');
end

% construct a weight vector based upon the quality of the initial fit of
% the ellipse to the pupil perimeter
if p.Results.nFrames ~= Inf
    if isfield(pupilData,p.Results.whichFitFieldError)
        errorWeights = (1./pupilData.(p.Results.whichFitFieldError)(p.Results.startFrame:p.Results.nFrames));
        errorWeights = errorWeights ./ nanmean(errorWeights);
    else
        errorWeights=ones(1,size(ellipses,1));
    end
else
    if isfield(pupilData,p.Results.whichFitFieldError)
        errorWeights = (1./pupilData.(p.Results.whichFitFieldError)(p.Results.startFrame:end));
        errorWeights = errorWeights ./ nanmean(errorWeights);
    else
        errorWeights=ones(1,size(ellipses,1));
    end
end

% define an anonymous function to measure mean error (the L1 norm)
errorFunc = @(x) nanmean( errorWeights.*ellipseCenterPredictionErrors(ellipses, x(1:3), x(4), p.Results.projectionModel) );

% define some search options
options = optimset('fmincon');
options = optimset(options,'Diagnostics','off','Display','off','LargeScale','off','Algorithm','interior-point');

% perform an initial fit
[bestFitSceneGeometry, fVal] = ...
    fmincon(errorFunc, x0, [], [], [], [], p.Results.sceneGeometryLowerBounds, p.Results.sceneGeometryUpperBounds, [], options);

% Split the ellipses by quantiles of distance error and then repeat the
% fit. This helps us estimate the range of acceptable scene parameters
distanceError = ellipseCenterPredictionErrors(ellipses, bestFitSceneGeometry(1:3), bestFitSceneGeometry(4), p.Results.projectionModel);
for nn = 1:p.Results.numberFitDivisions
    quantileLowBound = quantile(distanceError,(nn-1)/p.Results.numberFitDivisions);
    quantileUpperBound = quantile(distanceError,nn/p.Results.numberFitDivisions);
    withinBoundIdx = find((distanceError > quantileLowBound) .* (distanceError <= quantileUpperBound));
    errorFuncByQuantile = @(x) nanmean( errorWeights(withinBoundIdx).*ellipseCenterPredictionErrors(ellipses(withinBoundIdx,:), x(1:3), x(4), p.Results.projectionModel) );
    fitSceneGeometryByQuantile(nn,:) = ...
        fmincon(errorFuncByQuantile, bestFitSceneGeometry, [], [], [], [], p.Results.sceneGeometryLowerBounds, p.Results.sceneGeometryUpperBounds, [], options);    
end

% Derive the upper and lower bounds for the scene geometry, accounting for
% the possibility that the bestFitSceneGeometry has values that are higher
% or lower than any of the quantile split results.
lbSceneGeometrty = min([fitSceneGeometryByQuantile; bestFitSceneGeometry]);
ubSceneGeometrty = max([fitSceneGeometryByQuantile; bestFitSceneGeometry]);

% plot the results of the CoP estimation if requested
if ~isempty(p.Results.sceneDiagnosticPlotFileName)
    [~, predictedEllipseCenterXY] = ellipseCenterPredictionErrors(ellipses, bestFitSceneGeometry(1:3), bestFitSceneGeometry(4), p.Results.projectionModel);
    figHandle = figure('visible','off');
    plot(ellipses(:,1), ellipses(:,2), '.k')
    hold on
    plot(predictedEllipseCenterXY(:,1), predictedEllipseCenterXY(:,2), '.b')
    plot(x0(1),x0(2), 'xr')
    plot(bestFitSceneGeometry(1),bestFitSceneGeometry(2), 'og')
    xlim ([0 p.Results.sceneDiagnosticPlotSizeXY(1)])
    ylim ([0 p.Results.sceneDiagnosticPlotSizeXY(2)])
    set(gca,'Ydir','reverse')
    title('Estimate center of eye rotation from pupil ellipses')
    legend('ellipse centers','predicted ellipse centers', 'Most circular ellipse','Best fit CoR')
    saveas(figHandle,p.Results.sceneDiagnosticPlotFileName);
end

% assemble and save the sceneGeometry
sceneGeometry.eyeCenter.X = bestFitSceneGeometry(1);
sceneGeometry.eyeCenter.X_bounds = [lbSceneGeometrty(1) ubSceneGeometrty(1)];
sceneGeometry.eyeCenter.Y = bestFitSceneGeometry(2);
sceneGeometry.eyeCenter.Y_bounds = [lbSceneGeometrty(2) ubSceneGeometrty(2)];
sceneGeometry.eyeCenter.Z = bestFitSceneGeometry(3);
sceneGeometry.eyeCenter.Z_bounds = [lbSceneGeometrty(3) ubSceneGeometrty(3)];
sceneGeometry.eyeRadius = bestFitSceneGeometry(4);
sceneGeometry.eyeRadius_bounds = [lbSceneGeometrty(4) ubSceneGeometrty(4)];
sceneGeometry.meta = p.Results;
sceneGeometry.meta.units = 'pixelsOnTheScenePlane';
sceneGeometry.meta.meanDistanceError = fVal;

if ~isempty(sceneGeometryFileName)
    save(sceneGeometryFileName,'sceneGeometry');
end

% alert the user that we are done with the routine
if strcmp(p.Results.verbosity,'full')
    toc
    fprintf('\n');
end

end % main function




%% LOCAL FUCTIONS

function [distanceError, predictedEllipseCenterXY] = ellipseCenterPredictionErrors(ellipses, candidateEyeCenterOfRotation, eyeRadius, projectionModel)
% distanceError = ellipseCenterPredictionErrors(ellipses, candidateEyeCenterOfRotation, eyeRadius, projectionModel)
%
% Calculate the error in predicting the center of ellipses on the image
% plane based upon the eccentricity and theta of the ellipse and the
% candidate scene geometry
%
% input
%  ellipses - set of ellipses in transparent form.
%  candidateEyeCenterOfRotation -[X Y Z] coordinates of the center of the
%       eye
%  eyeRadius - radius of the eye in pixels
%

% initialize variables
distanceError = nan(size(ellipses,1),1);
predictedEllipseCenterXY = nan(size(ellipses,1),2);

% loop through ellipses
for ii = 1:size(ellipses,1)
    if ~any(isnan(ellipses(ii,:)))
        % Obtain the azimuth and elevation of the eye implied by the
        % eccentricity and theta of the ellipse, without any knowledge of
        % scene geometry. Two possible solutions will be returned. We nan
        % out the ellipse area, as this cannot inform the calculation
        pupilEllipseOnImagePlane = ellipses(ii,:);
        pupilEllipseOnImagePlane(3) = nan;
        [reconstructedPupilAzi, reconstructedPupilEle, ~] = pupilProjection_inv(pupilEllipseOnImagePlane, [nan nan nan], nan, projectionModel);
        
        % Obtain the x and y position of the projection of a pupil at this
        % azimuth and elevation onto the pupil plane, using the passed
        % candidate eye center. A nan value is passed for pupil area.
        projectedEllipse(1,:) = pupilProjection_fwd(reconstructedPupilAzi(1), reconstructedPupilEle(1), nan, candidateEyeCenterOfRotation, eyeRadius, projectionModel);
        projectedEllipse(2,:) = pupilProjection_fwd(reconstructedPupilAzi(2), reconstructedPupilEle(2), nan, candidateEyeCenterOfRotation, eyeRadius, projectionModel);
                
        % Calculate the distance between the projected ellipse center and
        % the observed ellipse center and save the smaller of the two
        % candidate calues
        euclideanDistances(1) = sqrt( (projectedEllipse(1,1) - pupilEllipseOnImagePlane(1))^2 + (projectedEllipse(1,2) - pupilEllipseOnImagePlane(2))^2 );
        euclideanDistances(2) = sqrt( (projectedEllipse(2,1) - pupilEllipseOnImagePlane(1))^2 + (projectedEllipse(2,2) - pupilEllipseOnImagePlane(2))^2 );
        [distanceError(ii), minIdx] = min(euclideanDistances);
        predictedEllipseCenterXY(ii,1:2) = projectedEllipse(minIdx,1:2);
    else
        distanceError(ii) = nan;
        predictedEllipseCenterXY(ii,1:2) = nan;
    end
end

end % local function
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
%	pupilFileName - full path to a pupilData file, or a cell array of such
%      paths
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
% Options (environment)
%   tbSnapshot - the passed tbSnapshot output that is to be saved along
%      with the data
%   timestamp / username / hostname - these are automatically derived and
%      saved within the p.Results structure.
%


%% input parser
p = inputParser; p.KeepUnmatched = true;

% required input
p.addRequired('pupilFileName',@(x)(iscell(x) | ischar(x)));
p.addRequired('sceneGeometryFileName',@ischar);

% Optional analysis params
p.addParameter('projectionModel','pseudoPerspective',@ischar);
p.addParameter('eyeRadius',125,@isnumeric);
p.addParameter('cameraDistanceInPixels',1200,@isnumeric);
p.addParameter('sceneGeometryLowerBounds',[0, 0, 1325, 25],@isnumeric);
p.addParameter('sceneGeometryUpperBounds',[640, 480, 1325, 250],@isnumeric);
p.addParameter('whichFitFieldMean','ellipseParamsUnconstrained_mean',@ischar);
p.addParameter('whichFitFieldError','ellipseParamsUnconstrained_rmse',@ischar);

% verbosity and plotting control
p.addParameter('verbosity', 'none', @isstr);
p.addParameter('sceneDiagnosticPlotFileName', [],@(x)(isempty(x) | ischar(x)));
p.addParameter('sceneDiagnosticPlotSizeXY', [640 480],@isnumeric);

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
if iscell(pupilFileName)
    ellipses = [];
    errorWeights = [];
    for cc = 1:length(pupilFileName)
        load(pupilFileName{cc})
        ellipses = [ellipses;pupilData.(p.Results.whichFitFieldMean)];
        errorWeights = [errorWeights; pupilData.(p.Results.whichFitFieldError)];
    end
else
    load(pupilFileName)
    ellipses = pupilData.(p.Results.whichFitFieldMean);
    errorWeights = pupilData.(p.Results.whichFitFieldError);
end

% find the most circular ellipse and use the X Y coordinate as the initial
% guess for the CoP XY coordinates
[~, minEccentricityIdx] = min(pupilData.(p.Results.whichFitFieldMean)(:,4));

% Set up the initial guess at the sceneGeometry.
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
errorWeights = 1./errorWeights;
errorWeights = errorWeights ./ nanmean(errorWeights);

% We divide the ellipse centers amongst a 2D set of bins. We will
% ultimately minimize the fitting error across bins
[ellipseCenterCounts,Xedges,Yedges,binXidx,binYidx] = ...
    histcounts2(ellipses(:,1),ellipses(:,2),10);

% anonymous functions for row and column identity given array position
rowIdx = @(b) fix( (b-1) ./ (size(ellipseCenterCounts,2)) ) +1;
colIdx = @(b) 1+mod(b-1,size(ellipseCenterCounts,2));

% Create a cell array of index positions corresponding to each of the 2D
% bins
idxByBinPosition = ...
    arrayfun(@(b) find( (binXidx==rowIdx(b)) .* (binYidx==colIdx(b)) ),1:1:numel(ellipseCenterCounts),'UniformOutput',false);

% anonymous function to return the weighted distance error for each ellipse
distanceErrorVector = @(x, b) errorWeights(b).*ellipseCenterPredictionErrors(ellipses(b,:), x(1:3), x(4), p.Results.projectionModel );

% anonymous function to return the median distance error in each bin
medianDistanceErrorByBin = @(x) ...
    cell2mat(cellfun(@(b) nanmedian(distanceErrorVector(x, b)), idxByBinPosition, 'UniformOutput', false));

% anonymous function to calculate RMSE error across bins
errorFunc = @(x) sqrt( nanmean(medianDistanceErrorByBin(x).^2) );


% define some search options
options = optimset('fmincon');
options = optimset(options,'Diagnostics','off','Display','off','LargeScale','off','Algorithm','interior-point');

% perform the fit
[bestFitSceneGeometry, fVal] = ...
    fmincon(errorFunc, x0, [], [], [], [], p.Results.sceneGeometryLowerBounds, p.Results.sceneGeometryUpperBounds, [], options);

% assemble and save the sceneGeometry
sceneGeometry.eyeCenter.X = bestFitSceneGeometry(1);
sceneGeometry.eyeCenter.Y = bestFitSceneGeometry(2);
sceneGeometry.eyeCenter.Z = bestFitSceneGeometry(3);
sceneGeometry.eyeRadius = bestFitSceneGeometry(4);
sceneGeometry.meta = p.Results;
sceneGeometry.meta.units = 'pixelsOnTheScenePlane';
sceneGeometry.meta.distanceError = fVal;

if ~isempty(sceneGeometryFileName)
    save(sceneGeometryFileName,'sceneGeometry');
end

% plot the results of the CoP estimation if requested
if ~isempty(p.Results.sceneDiagnosticPlotFileName)
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
    scatter(ellipses(:,1),ellipses(:,2),2,'filled', ...
       'MarkerFaceAlpha',1/8,'MarkerFaceColor',[0.5 0.5 0.5]);
    hold on
    % get the predicted ellipse centers
    [~, predictedEllipseCenterXY] = ellipseCenterPredictionErrors(ellipses, bestFitSceneGeometry(1:3), bestFitSceneGeometry(4), p.Results.projectionModel);
    % plot the predicted ellipse centers
    scatter(predictedEllipseCenterXY(:,1),predictedEllipseCenterXY(:,2),2,'filled', ...
       'MarkerFaceAlpha',2/8,'MarkerFaceColor',[0 0 1]);
    % plot the position of the most circular ellipse
    plot(x0(1),x0(2), 'xr', 'MarkerSize', 3);
    % plot the estimated center of rotation of the eye
    plot(bestFitSceneGeometry(1),bestFitSceneGeometry(2), 'og', 'MarkerSize', 3);
    % label and clean up the plot
    xlim ([Xedges(1)-2*binSpaceX Xedges(end)+2*binSpaceX]);
    ylim ([Yedges(1)-2*binSpaceY Yedges(end)+2*binSpaceY]);
    axis equal
    set(gca,'Ydir','reverse')
    title('Estimate center of eye rotation from pupil ellipses')
    
    % Create a legend
    hSub = subplot(2,2,2);
    scatter(nan, nan,1,'filled', ...
        'MarkerFaceAlpha',1/8,'MarkerFaceColor',[0.5 0.5 0.5]);
    hold on
    scatter(nan, nan,1,'filled', ...
        'MarkerFaceAlpha',1/8,'MarkerFaceColor',[0 0 1]);
    plot(nan, nan, 'xr', 'MarkerSize', 2);
    plot(nan, nan, 'og', 'MarkerSize', 2);
    set(hSub, 'Visible', 'off');
    legend({'ellipse centers','predicted ellipse centers', 'Most circular ellipse','Best fit CoR'});
    
    % Next, plot the ellipse counts and error values by bin
    subplot(2,2,3)    
    image = flipud(ellipseCenterCounts');
    image(image==0)=nan;
    [nr,nc] = size(image);
    pcolor([flipud(image) nan(nr,1); nan(1,nc+1)]);
    caxis([0 max(max(image))]);
    shading flat;
    axis equal
    % Set the axis backgroud to dark gray
    set(gcf,'Color',[1 1 1]); set(gca,'Color',[.75 .75 .75]); set(gcf,'InvertHardCopy','off');
    set(gca,'Ydir','reverse')
    c = colorbar;
    c.Label.String = 'Ellipse counts per bin';
    xticks(0.5:1:size(image,1)+.5);
    xticklabels(Xedges);
    yticks(0.5:1:size(image,2)+.5);
    yticklabels(Yedges);
    xlim([1 size(image,1)+1]);
    ylim([1 size(image,2)+1]);

    subplot(2,2,4)    
    image = flipud(reshape(medianDistanceErrorByBin(bestFitSceneGeometry),size(ellipseCenterCounts,2),size(ellipseCenterCounts,1)));
    [nr,nc] = size(image);
    pcolor([flipud(image) nan(nr,1); nan(1,nc+1)]);
    caxis([0 max(max(image))]);
    shading flat;
    axis equal
    % Set the axis backgroud to dark gray
    set(gcf,'Color',[1 1 1]); set(gca,'Color',[.75 .75 .75]); set(gcf,'InvertHardCopy','off');
    set(gca,'Ydir','reverse')
    c = colorbar;
    c.Label.String = 'Median distance error';
    xticks(0.5:1:size(image,1)+.5);
    xticklabels(Xedges);
    yticks(0.5:1:size(image,2)+.5);
    yticklabels(Yedges);
    xlim([1 size(image,1)+1]);
    ylim([1 size(image,2)+1]);
    
    saveas(figHandle,p.Results.sceneDiagnosticPlotFileName);
    close(figHandle)
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
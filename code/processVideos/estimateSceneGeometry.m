function initialSceneGeometry = estimateSceneGeometry(pupilFileName, sceneGeometryFileName, varargin)
% Estimate eye radius and eye center given a set of image plane ellipses
%
% Description:
%   This function searches over the set of ellipses in the passed pupil
%   file to estimate the sceneGeometry features in units of pixels on the
%   scene. The routine identifies the eye radius and the [X, Y, Z]
%   coordinates on the scene plane of the center of rotation of an eye. The
%   search attempts to minimize the error associated with the prediction of
%   the eccentricity and theta of ellipses in each of 100 "bins" of x, y
%   ellipse center position on the image plane.
%
%   Different projection models can be used to guide this calculation. The
%   orthogonal model assumes that ellipses  on the scene plane are
%   orthogonal projections of a circular pupil the center of which rotates
%   around the eye center. The pseudoPerspective model adjusts the x, y
%   position of the center of an ellipse on the image plane given the
%   increased distance from the image plane when the eye is rotated.
%
% 	Note: the search for both eye radius and eyeCenter.Z is not
% 	sufficiently constrainted. Therefore, the boundaries for one of these
% 	should be locked.
%
% Inputs:
%	pupilFileName         - Full path to a pupilData file, or a cell array
%                           of such paths.
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
%  'projectionModel'      - Options are: 'orthogonal' and 'perspective'
%  'sceneGeometryLB'      - Lower bounds for the sceneGeometry parameter
%                           search. This is a 4x1 vector specifying
%                           eyeCenter.X, eyeCenter.Y, eyeCenter.Z and
%                           eyeRadius.
%  'sceneGeometryUB'      - The corresponding upper bounds.
%  'cameraDistanceInPixels' - This is used (along with eyeRadius) to
%                           construct an initial guess for eyeCenter.Z
%  'eyeRadius             - Under the orthogonal projection case, this
%                           value is stored and used for all subsequent
%                           calculations. Under the pseudoPerspective case,
%                           this is the initial guess for eyeRadius.
%  'whichFitFieldMean'    - Identifies the field in pupilData that contains
%                           the ellipse fit params for which the search
%                           will be conducted.
%  'whichFitFieldError'   - Identifies the pupilData field that has error
%                           values for the ellipse fit params.
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

% Optional environment parameters
p.addParameter('tbSnapshot',[],@(x)(isempty(x) | isstruct(x)));
p.addParameter('timestamp',char(datetime('now')),@ischar);
p.addParameter('username',char(java.lang.System.getProperty('user.name')),@ischar);
p.addParameter('hostname',char(java.net.InetAddress.getLocalHost.getHostName),@ischar);

% Optional analysis params
p.addParameter('intrinsicCameraMatrix',[772.5483 0 320; 0 772.5483 240; 0 0 1],@isnumeric);
p.addParameter('extrinsicTranslationVector',[0; 0; 38],@isnumeric);
p.addParameter('extrinsicRotationMatrix',[1 0 0; 0 -1 0; 0 0 -1],@isnumeric);
p.addParameter('eyeRadius',12,@isnumeric);
p.addParameter('extrinsicTranslationVectorLB',[-5; -5; 20],@isnumeric);
p.addParameter('extrinsicTranslationVectorUB',[5; 5; 50],@isnumeric);
p.addParameter('eyeRadiusLB',11,@isnumeric);
p.addParameter('eyeRadiusUB',14,@isnumeric);
p.addParameter('ellipseConstraintTolerance',0.02,@isnumeric);

p.addParameter('whichEllipseFitField','initial',@ischar);

% parse
p.parse(pupilFileName, sceneGeometryFileName, varargin{:})


%% Announce we are starting
if strcmp(p.Results.verbosity,'full')
    tic
    fprintf(['Estimating scene geometry from pupil ellipses. Started ' char(datetime('now')) '\n']);
end


%% Identify the ellipses that will guide the sceneGeometry estimation
% load pupil data
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

% We divide the ellipse centers amongst a 5x5 set of bins across image
% space. We will ultimately minimize the fitting error across bins
[ellipseCenterCounts,Xedges,Yedges,binXidx,binYidx] = ...
    histcounts2(ellipses(:,1),ellipses(:,2),5);

% Anonymous functions for row and column identity given array position
rowIdx = @(b) fix( (b-1) ./ (size(ellipseCenterCounts,2)) ) +1;
colIdx = @(b) 1+mod(b-1,size(ellipseCenterCounts,2));

% Create a cell array of index positions corresponding to each of the 2D
% bins
idxByBinPosition = ...
    arrayfun(@(b) find( (binXidx==rowIdx(b)) .* (binYidx==colIdx(b)) ),1:1:numel(ellipseCenterCounts),'UniformOutput',false);

% Identify which bins are not empty
filledBinIdx = find(~cellfun(@isempty, idxByBinPosition));

% Identify the ellipses in each filled bin with the lowest fitting error
[~, idxMinErrorEllipseWithinBin] = arrayfun(@(x) min(ellipseFitSEM(idxByBinPosition{x})), filledBinIdx, 'UniformOutput', false);
returnTheMin = @(binContents, x)  binContents(idxMinErrorEllipseWithinBin{x});
ellipseArrayList = cellfun(@(x) returnTheMin(idxByBinPosition{filledBinIdx(x)},x),num2cell(1:1:length(filledBinIdx)));


%% Create the initial sceneGeometry structure and bounds
% sceneGeometry
initialSceneGeometry.eyeRadius = p.Results.eyeRadius;
initialSceneGeometry.intrinsicCameraMatrix = p.Results.intrinsicCameraMatrix;
initialSceneGeometry.extrinsicTranslationVector = p.Results.extrinsicTranslationVector;
initialSceneGeometry.extrinsicRotationMatrix = p.Results.extrinsicRotationMatrix;

% Initial search point
x0 = [initialSceneGeometry.extrinsicTranslationVector; initialSceneGeometry.eyeRadius];

% Bounds
lb = [p.Results.extrinsicTranslationVectorLB; p.Results.eyeRadiusLB];
ub = [p.Results.extrinsicTranslationVectorUB; p.Results.eyeRadiusUB];

%% Perform the search

% Call out to the local search as the nested function creates a static
% workspace, preventing us from loading the pupilData file in this
% function.
[sceneGeometry, rmseDistanceError, medianDistanceErrorPerBin] = ...
    performSearch(initialSceneGeometry, ellipses, ellipseArrayList, x0, lb, ub, p.Results.ellipseConstraintTolerance);

%% Save the sceneGeometry file

% add a meta field
sceneGeometry.meta = p.Results;
sceneGeometry.meta.rmseDistanceError = rmseDistanceError;

if ~isempty(sceneGeometryFileName)
    save(sceneGeometryFileName,'sceneGeometry');
end


%% Create a sceneGeometry plot
% plot the results of the CoP estimation if requested
if ~isempty(p.Results.sceneDiagnosticPlotFileName)
    
    figHandle = figure('visible','on');
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
    scatter(ellipses(ellipseArrayList,1),ellipses(ellipseArrayList,2),'o','filled', ...
       'MarkerFaceAlpha',2/8,'MarkerFaceColor',[0 0 0]);
    hold on

    % get the predicted ellipse centers
        [~, projectedEllipses] = ...
            arrayfun(@(x) pupilProjection_inv...
                    (...
                    ellipses(x,:),...
                    sceneGeometry,...
                    'constraintTolerance', p.Results.ellipseConstraintTolerance...
                    ),...
                ellipseArrayList,'UniformOutput',false);
    projectedEllipses=vertcat(projectedEllipses{:});

    % plot the projected ellipse centers
    scatter(projectedEllipses(:,1),projectedEllipses(:,2),'o','filled', ...
       'MarkerFaceAlpha',2/8,'MarkerFaceColor',[0 0 1]);

    % connect the centers with lines
    for ii=1:length(ellipseArrayList)
        plot([projectedEllipses(ii,1) ellipses(ellipseArrayList(ii),1)],[projectedEllipses(ii,2) ellipses(ellipseArrayList(ii),2)],'-r');
    end
    
   % plot the estimated center of rotation of the eye
    centerOfRotationEllipse = pupilProjection_fwd([0 0 2], sceneGeometry);
    plot(centerOfRotationEllipse(1),centerOfRotationEllipse(2), '+r', 'MarkerSize', 5);

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
    plot(nan, nan, '+r', 'MarkerSize', 2);
    set(hSub, 'Visible', 'off');
    legend({'ellipse centers','predicted ellipse centers', 'Best fit CoR'},'Location','southwestoutside');
    
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
    colorbar;
    title('Ellipse counts')
    xticks(1:1:size(image,1)+1);
    xticklabels(round(Xedges));
    xtickangle(90);
    yticks(1:1:size(image,2)+1);
    yticklabels(round(Yedges));
    xlim([1 size(image,1)+1]);
    ylim([1 size(image,2)+1]);

%     subplot(2,2,4)    
%     image = flipud(reshape(medianDistanceErrorByBin(bestFitSceneGeometry),size(ellipseCenterCounts,2),size(ellipseCenterCounts,1)));
%     [nr,nc] = size(image);
%     pcolor([flipud(image) nan(nr,1); nan(1,nc+1)]);
%     caxis([0 max(max(image))]);
%     shading flat;
%     axis equal
%     % Set the axis backgroud to dark gray
%     set(gcf,'Color',[1 1 1]); set(gca,'Color',[.75 .75 .75]); set(gcf,'InvertHardCopy','off');
%     set(gca,'Ydir','reverse')
%     colorbar;
%     title('Median distance error');
%     xticks(1:1:size(image,1)+1);
%     xticklabels(round(Xedges));
%     xtickangle(90);
%     yticks(1:1:size(image,2)+1);
%     yticklabels(round(Yedges));
%     xlim([1 size(image,1)+1]);
%     ylim([1 size(image,2)+1]);
    
    saveas(figHandle,p.Results.sceneDiagnosticPlotFileName);
    close(figHandle)
end

% alert the user that we are done with the routine
if strcmp(p.Results.verbosity,'full')
    toc
    fprintf('\n');
end

end % main function



function [adjustedSceneGeometry, rmseDistanceError, medianDistanceErrorByBin] = performSearch(sceneGeometry, ellipses, ellipseArrayList, x0, lb, ub, ellipseConstraintTolerance)

% Define search options
options = optimoptions(@patternsearch, ...
    'Display','iter',...
    'AccelerateMesh',false,...
    'FunctionTolerance',0.01);

% Define anonymous functions for the objective and constraint
objectiveFun = @objfun; % the objective function, nested below

% Define nested variables for within the search
medianDistanceErrorByBin=[];

[x, rmseDistanceError] = patternsearch(objectiveFun, x0,[],[],[],[],lb,ub,[],options);
    function fval = objfun(x)
        candidateSceneGeometry = sceneGeometry;
        candidateSceneGeometry.extrinsicTranslationVector = x(1:3);
        candidateSceneGeometry.eyeRadius = x(4);
        [~, ~, centerErrors] = ...
            arrayfun(@(x) pupilProjection_inv...
                    (...
                    ellipses(x,:),...
                    candidateSceneGeometry,...
                    'constraintTolerance', ellipseConstraintTolerance...
                    ),...
                ellipseArrayList,'UniformOutput',false);
        
        % Now compute objective function as the RMSE of the distance
        % between the taget and modeled ellipses
        medianDistanceErrorByBin = cellfun(@median, centerErrors);
        fval = sqrt(mean(medianDistanceErrorByBin.^2));
    end

adjustedSceneGeometry = sceneGeometry;
adjustedSceneGeometry.extrinsicTranslationVector = x(1:3);
adjustedSceneGeometry.eyeRadius = x(4);

end % local search function



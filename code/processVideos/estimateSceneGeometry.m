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
%   eyeRadiusInPixels - an assigned valye for eye radius. When operating 
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
p.addParameter('projectionModel','orthogonal', @ischar);
p.addParameter('eyeRadiusInPixels',250,@isnumeric);
p.addParameter('CoRLowerBound',[-500, -500, 0],@isnumeric);
p.addParameter('CoRUpperBound',[1000, 1000, 1000],@isnumeric);
p.addParameter('whichFitField','orthogonal', @ischar);
p.addParameter('whichFitFieldMean','ellipseParamsUnconstrained_mean',@ischar);
p.addParameter('whichFitFieldError','ellipseParamsUnconstrained_rmse',@ischar);

% verbosity and plotting control
p.addParameter('verbosity', 'none', @isstr);
p.addParameter('sceneDiagnosticPlotFileName', [],@(x)(isempty(x) | ischar(x)));

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
            0];
    case 'perspective'
        error('not yet implemented');
end

% construct a weight vector based upon the quality of the initial fit of
% the ellipse to the pupil perimeter
if p.Results.nFrames ~= Inf
    if isfield(pupilData,p.Results.whichFitFieldError)
        errorWeights = (1./pupilData.(p.Results.whichFitFieldError)(p.Results.startFrame:p.Results.nFrames));
    else
        errorWeights=ones(1,size(ellipses,1));
    end
else
    if isfield(pupilData,p.Results.whichFitFieldError)
        errorWeights = (1./pupilData.(p.Results.whichFitFieldError)(p.Results.startFrame:end));
    else
        errorWeights=ones(1,size(ellipses,1));
    end
end

% define an anonymous function to measure root mean squared error (RMSE)
errorFunc = @(x) sqrt(nanmean(errorWeights.*distanceToCandidateEyeCenterOfRotation(ellipses, x, p.Results.eyeRadiusInPixels, 'projectionModel', p.Results.projectionModel).^2));

% define some search options
options = optimset('fmincon');
options = optimset(options,'Diagnostics','off','Display','off','LargeScale','off','Algorithm','interior-point');

% perform the fit
[bestFitCoR, fVal] = fmincon(errorFunc, x0, [], [], [], [], p.Results.CoRLowerBound, p.Results.CoRUpperBound, [], options);

% plot the results of the CoP estimation if requested
if ~isempty(p.Results.sceneDiagnosticPlotFileName)
    [~, ellipsesCoRs] = distanceToCandidateEyeCenterOfRotation(ellipses,bestFitCoR,p.Results.eyeRadiusInPixels,'projectionModel', p.Results.projectionModel);
    figHandle = figure('visible','off');
    plot(ellipses(:,1), ellipses(:,2), '.k')
    hold on
    plot(ellipsesCoRs(:,1), ellipsesCoRs(:,2), '.b')
    plot(x0(1),x0(2), 'xr')
    plot(bestFitCoR(1),bestFitCoR(2), 'og')
    title('Estimate Center of Rotation from pupil ellipses')
    legend('ellipse centers','CoR from each ellipse', 'Most circular ellipse','Best fit CoR')
    saveas(figHandle,p.Results.sceneDiagnosticPlotFileName);
end

% assemble and save the sceneGeometry
sceneGeometry.eyeCenter.X = bestFitCoR(1);
sceneGeometry.eyeCenter.Y = bestFitCoR(2);
sceneGeometry.eyeCenter.Z = bestFitCoR(3);
sceneGeometry.eyeCenter.RMSE = fVal;
sceneGeometry.eyeRadius = p.Results.eyeRadiusInPixels;
sceneGeometry.meta = p.Results;
sceneGeometry.meta.units = 'pixelsOnTheScenePlane';

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

function [distances, ellipsesCoRs] = distanceToCandidateEyeCenterOfRotation(ellipses, candidateEyeCenterOfRotation, eyeRadiusInPixels, varargin)
%  [distances, ellipsesCoRs] = distanceToCandidateEyeCenterOfRotation(ellipses, candidateEyeCenterOfRotation, eyeRadiusInPixels, varargin)
%
% this function finds the distance of candidate center of rotation for
% each ellipses and the center of projection on the scene. Ellipses and
% candidate eye center of rotation must be in the same units.
% (e.g. Pixels)
%
% input
%  ellipses - set of ellipses in transparent form.
%  candidateEyeCenterOfRotation -[X Y Z] coordinates of the center. If
%       orthogonal projection, the Z coordinate (distance from the scene plane,
%       will be set to zero).
%  eyeRadiusInPixels - radius of the eye in pixels
%
%  Output
%     distances
%     ellipsesCoRs

%% input parser
p = inputParser; p.KeepUnmatched = true;

% required input
p.addRequired('ellipses',@isnumeric);
p.addRequired('candidateEyeCenterOfRotation',@isnumeric);
p.addRequired('eyeRadiusInPixels',@isnumeric);

% Analysis parameters
p.addParameter('projectionModel','orthogonal',@ischar);

% parse
p.parse(ellipses, candidateEyeCenterOfRotation, eyeRadiusInPixels, varargin{:})

%% main

% initialize variables
distances = nan(size(ellipses,1),2);
ellipsesCoRs = nan(size(ellipses,1),3);

switch p.Results.projectionModel
    case 'orthogonal'
        
        % loop through ellipses, find all distances and pick center of projection
        for ii = 1:size(ellipses,1)
            if ~any(isnan(ellipses(ii,:)))
                % find candidates Center of rotation
                ellipseCoRCandidates = calcEyeCenterOfRotation(ellipses(ii,:), eyeRadiusInPixels, 'projectionModel', p.Results.projectionModel);
                
                % find Eucledian distance for each ellipse CoR
                for jj = 1:size(ellipseCoRCandidates,1)
                    distances(ii,jj) = ...
                        sqrt((candidateEyeCenterOfRotation(1) - ellipseCoRCandidates(jj,1))^2 + (candidateEyeCenterOfRotation(2) - ellipseCoRCandidates(jj,2))^2 + (candidateEyeCenterOfRotation(3) - ellipseCoRCandidates(jj,3))^2);
                end
                % select ellipseCoR with the min distance from the scene CoR and store
                % it
                [~,minDistIDX] =min(distances(ii,:));
                ellipsesCoRs(ii,:) = ellipseCoRCandidates(minDistIDX,:);
            else
                continue
            end
        end
        
        distances = min(distances,[],2);
        
    case 'perspective'
        error('not implemented yet');
end % switch

end % local function
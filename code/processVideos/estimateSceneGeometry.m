function sceneGeometry = estimateSceneGeometry(pupilFileName, sceneGeometryFileName, varargin)
% sceneGeometry = estimateSceneGeometry(pupilFileName, sceneGeometryFileName)
%
% This function performs a search over the available ellipses to determine
% the sceneGeometry features in units of pixels on the scene. The returned
% features are the following:
%       eyeballCenter.X : X coordinate of the eyeball center (i.e. the assumed
%           center of rotation of the pupil) on the scene plane.
%       eyeballCenter.Y : Y coordinate of the eyeball center (i.e. the assumed
%           center of rotation of the pupil) on the scene plane.
%       eyeballCenter.Z : the orthogonal distance for the eyeball center from the
%           scene plane.
%       eyeballRadius : radius of the eyeball in pixels
%
% The routine can assume that the ellipses are orthogonal projections onto the
% scene plane of a circle whose center rotates around the eyeball center
% onto the scene plane.
% The most circular ellipse that is available provides the initial guess
% for the X Y coordinates of the center of the eyeball. 
% For each ellipse, the 2 candidate eyeballs are calculated.
% The routine will then look for the [X Y Z] coordinates that minimize the
% distance between the all the candidate center of projections.
%
% Note that at this stage, the best fit for the Z coordinate is not
% phisically sound. This might be due to the simplified assumptions of
% orthogonal projection, eyeball sphericity and absence of noise.
%
% Output
%  sceneGeometry - a structure with the fields eyeCenterX, Y, and Z, and
%  the field meta. The eyeCenter coordinates are in units of pixels on the
%  image plane.
%
% Input
%   pupilFileName -
%   sceneGeometryFileName -  name of the sceneGeometry file in which to
%       save the output.
%
%% input parser
p = inputParser; p.KeepUnmatched = true;

% required input
p.addRequired('pupilFileName',@isstr);
p.addRequired('sceneGeometryFileName',@isstr);

% optional inputs
p.addParameter('orthogonalProjection',true, @islogical);
p.addParameter('eyeballRadiusInPixels',250,@isnumeric);

% verbosity and plotting control
p.addParameter('verbosity', 'none', @isstr);
p.addParameter('plotResults', true, @islogical);

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

%% load pupil data

load(pupilFileName)
if p.Results.nFrames ~= Inf
    ellipses = pupilData.pInitialFitTransparent(p.Results.startFrame:p.Results.nFrames,:);
else
    ellipses = pupilData.pInitialFitTransparent(p.Results.startFrame:end,:);
end

%% find the most circular ellipse and use the X Y coordinate as the initial guess for the CoP XY coordinates

[minEccentricity, minEccentricityIDX] = min(pupilData.pInitialFitTransparent(:,4));

% under orthogonal projection hypotesis, the first guess for the CoP is the [X Y] coordinates of the most circular ellipse on
% the scene, the initial eyeball radius guess is user defined, and the Z coordinate is zero 
if p.Results.orthogonalProjection
    x0 = [pupilData.pInitialFitTransparent(minEccentricityIDX,1) ...
        pupilData.pInitialFitTransparent(minEccentricityIDX,2) ...
        0];
end
% add weights to ellipses based on how good is the fit
if p.Results.nFrames ~= Inf
    if isfield(pupilData,'fitError')
        errorWeights = (1./pupilData.fitError(p.Results.startFrame:p.Results.nFrames))';
    else
        errorWeights=ones(1,size(ellipses,1));
    end
else
    if isfield(pupilData,'fitError')
        errorWeights = (1./pupilData.fitError(p.Results.startFrame:end))';
    else
        errorWeights=ones(1,size(ellipses,1));
    end
end

% define an anonymous function to measure SSQ error
errorFunc = @(x) sqrt(nansum(errorWeights.*distanceToCandidateEyeBallCenter(ellipses, x, p.Results.eyeballRadiusInPixels, 'orthogonalProjection', p.Results.orthogonalProjection).^2));
bestFitCoP = fmincon(errorFunc, x0);


%% if so requested, plot the results of the CoP estimation

if p.Results.plotResults
    [distances, ellipsesCOPs] = distanceToCandidateEyeBallCenter(ellipses,bestFitCoP,p.Results.eyeballRadiusInPixels,'orthogonalProjection', p.Results.orthogonalProjection);
    figure
    plot(ellipsesCOPs(:,1), ellipsesCOPs(:,2), '.')
    hold on
    plot(x0(1),x0(2), 'xr')
    plot(bestFitCoP(1),bestFitCoP(2), 'og')
    title('Estimate Position of the Center of Projeciton for fitted ellipses')
    legend('Ellipses COP', 'COP of the most circular ellipse (reference COP)','Best fit CoP')
    xlim([0 320] * 2)
    ylim([0 240] * 2)
end

%% assemble the sceneGeometry and save it out
sceneGeometry.eyeballCenter.X = bestFitCoP(1);
sceneGeometry.eyeballCenter.Y = bestFitCoP(2);
sceneGeometry.eyeballCenter.Z = bestFitCoP(3); % meant as the distance from the scene plane
sceneGeometry.eyeballRadius = p.Results.eyeballRadiusInPixels;
sceneGeometry.units = 'pixelsOnTheScenePlane';
sceneGeometry.meta = p.Results;

if ~isempty(sceneGeometryFileName)
    save(sceneGeometryFileName,'sceneGeometry');
end

end % function


function [distances, ellipsesCoRs] = distanceToCandidateEyeBallCenter(ellipses, candidateEyeballCenter, eyeballRadiusInPixels, varargin)
%  [distances,ellipsesCoPs] = distanceToCandidateEyeBallCenter(ellipses,candidateEyeball,varargin)
%
% this function finds the distance of candidate center of rotation for
% each ellipses and the center of projection on the scene. Ellipses and
% candidate Eyeball must be in the same units.
% (e.g. Pixels)
%
% input 
%  ellipses - set of ellipses in transparent form.
%  candidateEyeballCenter -[X Y Z] coordinates of the center. If
%       orthogonal projection, the Z coordinate (distance from the scene plane,
%       will be set to zero).
%  eyeballRadiusInPixels - radius of the eyeball in pixel
% 
%  Output
%     distances
%     ellipsesCoRs
%% input parser
p = inputParser; p.KeepUnmatched = true;

% required input
p.addRequired('ellipses',@isnumeric);
p.addRequired('candidateEyeballCenter',@isnumeric);
p.addRequired('eyeballRadiusInPixels',@isnumeric);

% Analysis parameters
p.addParameter('orthogonalProjection',true,@islogical);

% parse
p.parse(ellipses, candidateEyeballCenter, eyeballRadiusInPixels, varargin{:})

if p.Results.orthogonalProjection
    
%% initialize variables
distances = nan(size(ellipses,1),2);
ellipsesCoRs = nan(size(ellipses,1),3);

%% loop through ellipses, find all distances and pick center of projection
for ii = 1:size(ellipses,1)
    if ~any(isnan(ellipses(ii,:)))
        % find candidates Center of rotation
        ellipseCoRCandidates = findCandidatesEyeball(ellipses(ii,:), eyeballRadiusInPixels, 'orthogonalProjection', p.Results.orthogonalProjection);
        
        % find eucledian distance for each ellipse COP
        for jj = 1:size(ellipseCoRCandidates,1)
            distances(ii,jj) = sqrt((candidateEyeballCenter(1) - ellipseCoRCandidates(jj,1))^2 + (candidateEyeballCenter(2) - ellipseCoRCandidates(jj,2))^2 + (candidateEyeballCenter(3) - ellipseCoRCandidates(jj,3))^2);
        end
        % select ellipseCoR with the min distance from the scene CoR and store
        % it
        [~,minDistIDX] =min(distances(ii,:));
        ellipsesCoRs(ii,:) = ellipseCoRCandidates(minDistIDX,:);
    else
        continue
    end
end

distances = min(distances');

else
    % Dev PLACEHOLDER - develop perspective correction case
end
end %function
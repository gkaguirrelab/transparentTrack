function sceneGeometry = estimateSceneGeometry(pupilFileName, sceneGeometryFileName, varargin)
% sceneGeometry = estimateSceneGeometry(pupilFileName, sceneGeometryFileName)
%
% This function performs a search over the available ellipses to determine
% the sceneGeometry features in units of pixels on the scene. The returned
% features are the following:
%       eyeball.X : X coordinate of the eyeball center (i.e. the assumed
%       center of rotation of the pupil) on the scene plane.
%       eyeball.Y : Y coordinate of the eyeball center (i.e. the assumed
%       center of rotation of the pupil) on the scene plane.
%       eyeball.Z : the orthogonal distance for the eyeball center from the
%       scene plane.
% 
% The routine assumes that the ellipses are orthogonal projections onto the
% scene plane of a circle whose center rotates around the eyeball center
% onto the scene plane.
% The most circular ellipse that is available provides the initial guess
% for the X Y coordinates of the center of the eyeball. The initial guess
% for the eyeball Z coordinate is arbitrary and can be customized.
% For each ellipse, the 2 candidate center of projection are calculated.
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
p.addParameter('initialEyeballRadiusPX',1000, @isnumeric);
p.addParameter('plotResults', true, @islogical);
p.addParameter('eyeRadiusOfCurvatureMm', 12, @islogical);

% verbosity
p.addParameter('verbosity', 'none', @isstr);

% flow control
p.addParameter('nFrames',Inf,@isnumeric);
p.addParameter('startFrame',1,@isnumeric);

%parse
p.parse(pupilFileName, sceneGeometryFileName, varargin{:})

%% load pupil data

load(pupilFileName)
if p.Results.nFrames ~= Inf
    ellipses = pupilData.pPosteriorMeanTransparent(p.Results.startFrame:p.Results.nFrames,:);
else
    ellipses = pupilData.pPosteriorMeanTransparent(p.Results.startFrame:end,:);
end

%% find the most circular ellipse and use the X Y coordinate as the initial guess for the CoP XY coordinates

[minEccentricity, minEccentricityIDX] = min(pupilData.pPosteriorMeanTransparent(:,4));

% first guess for the CoP is the location of the most circular ellipse on
% the scene and the initial eyeball radius guess as a Z coordinate
x0 = [pupilData.pPosteriorMeanTransparent(minEccentricityIDX,1) ...
    pupilData.pPosteriorMeanTransparent(minEccentricityIDX,2) ...
    p.Results.initialEyeballRadiusPX];

% add weights to ellipses based on how good is the fit
if p.Results.nFrames ~= Inf
errorWeights = (1./pupilData.fitError(p.Results.startFrame:p.Results.nFrames))';
else
    errorWeights = (1./pupilData.fitError(p.Results.startFrame:end))';
end

% define an anonymous function to measure SSQ error
errorFunc = @(x) sqrt(nansum(errorWeights.*distanceToCandidateEyeBallCenter(ellipses,x).^2));
bestFitCoP = fmincon(errorFunc, x0);


%% if so requested, plot the results of the CoP estimation

if p.Results.plotResults
    [distances, ellipsesCOPs] = distanceToCandidateEyeBallCenter(ellipses,bestFitCoP);
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
% note that the "eyeball" here is really the center of rotation
sceneGeometry.eyeball.X = bestFitCoP(1);
sceneGeometry.eyeball.Y = bestFitCoP(2);
sceneGeometry.eyeball.Z = bestFitCoP(3); % meant as the distance from the scene plane
sceneGeometry.units = 'pixelsOnTheScenePlane';
sceneGeometry.meta = p.Results;

save(sceneGeometryFileName,'sceneGeometry');
end % function
function sceneGeometry = estimateSceneGeometry(pupilFileName, sceneGeometryFileName, varargin)
% sceneGeometry = estimateSceneGeometry(pupilFileName, sceneGeometryFileName, varargin)
%
% HEADER
% Output
%  sceneGeometry - a structure with the fields eyeCenterX, Y, and Z, and
%  the field meta. The eyeCenter coordinates are in units of pixels on the
%  image plane.

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

% assemble the sceneGeometry structure here


end % function
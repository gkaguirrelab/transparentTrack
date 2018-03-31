function [ eyePose, bestMatchEllipseOnImagePlane, centerError, shapeError, areaError ] = inverseSearchWrapper( pupilEllipseOnImagePlane, sceneGeometry, eyePoseLB, eyePoseUB )
% Multistart search wrapper for pupilProjection_inv to avoid local minima
%
% Syntax:
%  [ eyePose, bestMatchEllipseOnImagePlane, centerError, shapeError, areaError ] = inverseSearchWrapper( pupilEllipseOnImagePlane, sceneGeometry, eyePoseLB, eyePoseUB )
%
% Description:
%   The pupilProjection_inv function is a search across forward model
%   solutions. It can find local minima. This wrapper performs multiple 
%   pupilProjection_inv searches, retaining the best solution.
%
% Inputs:
%   pupilEllipseOnImagePlane - A 1x5 vector that contains the parameters of
%                           pupil ellipse on the image plane cast in
%                           transparent form
%   sceneGeometry         - Structure. SEE: createSceneGeometry
%   eyePoseLB/UB          - A 1x4 vector that provides the lower (upper)
%                           bounds on the eyePose [azimuth, elevation,
%                           torsion, pupil radius].
%
% Outputs:
%   SEE: pupilProjection_inv for definitions
%


%% Set up variables

% This is out stopping criterion for a good inverse solution
centerErrorThreshold=1e-4;

% The first search will pass x0 as empty, allowing pupilProjection_inv to
% take its best guess as to the inverse solution. If stopping criteria are
% not met, then successive searches apply the azi and ele shift values to
% the eyePose returned by the initial search. The last of the azi and ele
% shift values are infinite, which will prompt pupilProjection_inv to
% select a random x0 starting point from within the eyePose bounds.
x0 = [];
aziShift = [0,-1,1,-1,1,Inf];
eleShift = [0,-1,1,1,-1,Inf];

% Set up a counter and a place to hold the best solution
searchCount = 0;
centerError = 1e6;


%% Perform the multi-start search
while 1
    searchCount = searchCount+1;
    % If we are not on the first search, use the initial search solution as
    % a starting point, adjusted by the values in azi and ele shift
    if searchCount > 1
        x0 = eyePoseInitial;
        x0(1) = x0(1) + aziShift(searchCount);
        x0(2) = x0(2) + eleShift(searchCount);
    end
    % Call out to pupilProjection_inv
	[eyePoseR, bestMatchEllipseOnImagePlaneR, centerErrorR, shapeErrorR, areaErrorR] = ...
        pupilProjection_inv(pupilEllipseOnImagePlane, sceneGeometry,'eyePoseLB',eyePoseLB,'eyePoseUB',eyePoseUB,'x0',x0);
    % If we just completed the first search, store this initial eyePose
    % result
    if searchCount == 1
        eyePoseInitial = eyePoseR;
    end
    % Did this search yield the best centerError value we have seen thus
    % far? If so, store outputs of the search so that we can return them
    if centerErrorR < centerError
        eyePose = eyePoseR;
        bestMatchEllipseOnImagePlane = bestMatchEllipseOnImagePlaneR;
        centerError = centerErrorR;
        shapeError = shapeErrorR;
        areaError = areaErrorR;
    end
    % Test if we are done our multi-start search. Our stopping criteria are
    % if we are below criterion for centerError or we have run out of
    % shifts to try for the starting position.
    if centerError<centerErrorThreshold || searchCount==length(aziShift)
        break
    end
end

end
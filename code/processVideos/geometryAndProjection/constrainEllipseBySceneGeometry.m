function [c, ceq]=constrainEllipseBySceneGeometry(transparentEllipseParams, sceneGeometry, constraintMarginEccenMultiplier, constraintMarginThetaDegrees, varargin)
% [c, ceq]=constrainEllipseBySceneGeometry(transparentEllipseParams, sceneGeometry, constraintMarginsEccenTheta, varargin)
%
% This function implements a non-linear constraint upon the ellipse fit
% to the pupil boundary. The goal of the limit is to constrain the
% eccentricity and theta of the ellipse to be a close match to that
% predicted by the x, y location of the center of the ellipse, given the
% sceneGeometry parameters.


% This function returns the expected eccentricity and tilt for an ellipse,
% given the location of the center and the scene geometry. The function can
% either use the Z coordinate for the eye in the scene geometry as the
% distance of the eye from the scene, or a user input distance range in
% pixel (2 elements vector). In the latter case, the routine will return
% the range of expected eccentricity and a single theta.
%
% Input (required)
% ellipseCenter - [X Y] coordinate for the ellipse center. Note that
%   passing the full parametrization of a transparent ellipse will work as
%   well.
% sceneGeometry - struct with scene geometry
%
% Output
% eccentricity - either a one element or 2 element vector with the expected
% 	eccentricity value or range.
% theta -the expected tilt value for the ellipse (does not depend on the
%   distance).
%


%% input parser
p = inputParser;

% required input
p.addRequired('transparentEllipseParams',@isnumeric);
p.addRequired('sceneGeometry',@isstruct);
p.addRequired('constraintMarginEccenMultiplier',@isnumeric);
p.addRequired('constraintMarginThetaDegrees',@isnumeric);

% optional analysis params
p.addParameter('projectionModel','orthogonal', @ischar);

%parse
p.parse(transparentEllipseParams, sceneGeometry, constraintMarginEccenMultiplier, constraintMarginThetaDegrees, varargin{:})


%% Calculate the predicted eccentricty and theta for this pupil center


% derive pupil center cartesian 3D coordinates
pupilCenter.Z = ...
    sqrt(...
        sceneGeometry.eyeRadius^2 - ...
        (transparentEllipseParams(1) - sceneGeometry.eyeCenter.X)^2 - ...
        (transparentEllipseParams(2) - sceneGeometry.eyeCenter.Y)^2 ...
        ) + ...
    sceneGeometry.eyeCenter.Z;
pupilAzi = atand((transparentEllipseParams(1) - sceneGeometry.eyeCenter.X)/(pupilCenter.Z - sceneGeometry.eyeCenter.Z));
pupilEle = atand((transparentEllipseParams(2) - sceneGeometry.eyeCenter.Y)/(pupilCenter.Z - sceneGeometry.eyeCenter.Z));

predictedTransparentEllipse = ...
    pupilProjection_fwd(pupilAzi, pupilEle, ...
    [sceneGeometry.eyeCenter.X sceneGeometry.eyeCenter.Y sceneGeometry.eyeCenter.Z],...
    'projectionModel',p.Results.projectionModel);

% calculate the range of allowed eccentricty and theta values
minAllowedEccentricity = predictedTransparentEllipse(4) ./ p.Results.constraintMarginEccenMultiplier;
maxAllowedEccentricity = predictedTransparentEllipse(4) .* p.Results.constraintMarginEccenMultiplier;

minAllowedTheta = predictedTransparentEllipse(5) - p.Results.constraintMarginThetaDegrees;
maxAllowedTheta = predictedTransparentEllipse(5) + p.Results.constraintMarginThetaDegrees;

% First constraint
%  Set ceq to zero only if the eccentricty of the ellipse to be tested is
%  within min / max allowed range
ceq = double(...
            ~((transparentEllipseParams(4) > minAllowedEccentricity) .* ...
            (transparentEllipseParams(4) < maxAllowedEccentricity))...
            );
% Second constraint
%  Set cq to zero only if the theta of the ellipse to be tested is
%  within min / max allowed range
c  = double(...
            ~((transparentEllipseParams(5) > minAllowedTheta) .* ...
            (transparentEllipseParams(5) < maxAllowedTheta))...
            );

end

function [c, ceq, projectedEllipseOnImagePlane]=constrainEllipseBySceneGeometry(pupilEllipseOnImagePlane, sceneGeometry, constraintFactor)
% [c, ceq]=constrainEllipseBySceneGeometry(transparentEllipseParams, sceneGeometry)
%
% This function implements a non-linear constraint upon the ellipse fit to
% the pupil boundary. The goal of the limit is to constrain the
% eccentricity and theta of the ellipse to be a close match to that
% predicted by the x, y location of the center of the ellipse, given the
% sceneGeometry parameters.
%
% Input (required)
%   pupilEllipseOnImagePlane - parameters of the ellipse to be fit cast in
%       transparent form
%   sceneGeometry - struct with scene geometry
%   constraintFactor - a positive multiplicative factor that adjusts how 
%       stingently the nonlinear constraint is evaluated. Values > 1
%       represent greater stringency, while values < 1 produce reduced
%       stringency.
%
% Output
%   c, ceq - the values of the nonlinear constraint function
%


% Extract some values from the sceneGeometry file
eyeCenterOfRotation = [sceneGeometry.eyeCenter.X sceneGeometry.eyeCenter.Y sceneGeometry.eyeCenter.Z];
eyeRadius = sceneGeometry.eyeRadius;
projectionModel = sceneGeometry.meta.projectionModel;

% Calculate the predicted x and y center for this eccentricity and theta
pupilEllipseOnImagePlane(3) = nan;
[reconstructedPupilAzi, reconstructedPupilEle, ~] = pupilProjection_inv(pupilEllipseOnImagePlane, eyeCenterOfRotation, eyeRadius, projectionModel);

% Obtain the x and y position of the projection of a pupil at this
% azimuth and elevation onto the pupil plane, using the passed
% sceneGeometry.
projectedEllipseOnImagePlane = pupilProjection_fwd(reconstructedPupilAzi, reconstructedPupilEle, nan, eyeCenterOfRotation, eyeRadius, projectionModel);

% Break out the predicted X and Y locations
predictedX = projectedEllipseOnImagePlane(1);
predictedY = projectedEllipseOnImagePlane(2);

% First constraint
%  Ceq reflects the Euclidean distance between the predicted and passed
%  center of the ellipse
c = constraintFactor .* sqrt( (pupilEllipseOnImagePlane(1) - predictedX).^2 + (pupilEllipseOnImagePlane(2) - predictedY).^2  );

% Second constraint
%  unused
ceq = [];

end

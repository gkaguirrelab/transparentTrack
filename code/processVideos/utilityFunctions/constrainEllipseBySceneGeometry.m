function [c, ceq]=constrainEllipseBySceneGeometry(pupilEllipseOnImagePlane, sceneGeometry)
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
%
% Output
%   c, ceq - the values of the nonlinear constraint function
%


% This is the value that is given to a violation of the constraint. This
% number needs to very large so that the constraint is respected despite
% large values in the objective function
constraintFactor = 1e20;

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

% First constraint
%  Ceq reflects the Euclidean distance between the predicted and passed
%  center of the ellipse
ceq = constraintFactor .* sqrt( (pupilEllipseOnImagePlane(1) - projectedEllipseOnImagePlane(1)).^2 + (pupilEllipseOnImagePlane(2) - projectedEllipseOnImagePlane(2)).^2  );

% Second constraint
%  unused
c = [];

end

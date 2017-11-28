function [c, ceq, projectedEllipseOnImagePlane]=constrainEllipseBySceneGeometry(pupilEllipseOnImagePlane, sceneGeometry)
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

% We increase the weight of the constraint by this factor to force
% compliance with the constraint in fitting
constraintFactor = 1e6;

% Extract some values from the sceneGeometry file
eyeCenterOfRotation = [sceneGeometry.eyeCenter.X sceneGeometry.eyeCenter.Y sceneGeometry.eyeCenter.Z];
eyeRadius = sceneGeometry.eyeRadius;
projectionModel = sceneGeometry.meta.projectionModel;

% Calculate the predicted x and y center for this eccentricity and theta
[reconstructedPupilAzi, reconstructedPupilEle, reconstructedPupilArea] = ...
    pupilProjection_inv(pupilEllipseOnImagePlane, eyeCenterOfRotation, eyeRadius, projectionModel);

% Calculate the projected ellipse on the image plane for the reconstructed
% azimuth and elevation
projectedEllipseOnImagePlane = ...
    pupilProjection_fwd(reconstructedPupilAzi, reconstructedPupilEle, reconstructedPupilArea, eyeCenterOfRotation, eyeRadius, projectionModel);

% First constraint
%  Ceq reflects the Euclidean distance between the predicted and passed
%  center of the ellipse
c = constraintFactor * sqrt( (pupilEllipseOnImagePlane(1) - projectedEllipseOnImagePlane(1)).^2 + (pupilEllipseOnImagePlane(2) - projectedEllipseOnImagePlane(2)).^2  );

% Second constraint
%  unused
ceq = [];

end

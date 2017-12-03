function [c, ceq, projectedEllipseOnImagePlane]=constrainEllipseBySceneGeometry(pupilEllipseOnImagePlane, sceneGeometry)
% A non-linear ellipse fit constraint function based on scene geometry
%
% Description:
%	This function implements a non-linear constraint upon an ellipse fit
%	performed within fmincon. The goal of the limit is to constrain the
%   eccentricity and theta of the ellipse to be a close match to that
%   predicted by the x, y location of the center of the ellipse, given the
%   sceneGeometry parameters.
%
% Input:
%   pupilEllipseOnImagePlane  - Parameters of an ellipse in transparent form
%                               to be evaluated
%   sceneGeometry             - Struct with scene geometry
%
% Output:
%   c     - The Euclidean distance between the X, Y center of the ellipse
%           and the center that would be predicted by the eccentricity and
%           theta of the ellipse given the scene constraints
%   ceq   - Unusued, so set to empty
%   projectedEllipseOnImagePlane - The transparent parameters of the
%           ellipse that is implied by the eccentricty and theta of the
%           input ellipse. While not used in the fitting routine, this
%           output is used by other routines to show the quality of the fit
%


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

% First constraint - Ceq reflects the Euclidean distance between the
%   predicted and passed center of the ellipse
c = sqrt( (pupilEllipseOnImagePlane(1) - projectedEllipseOnImagePlane(1)).^2 + (pupilEllipseOnImagePlane(2) - projectedEllipseOnImagePlane(2)).^2  );

% Second constraint - unused
ceq = [];

end % function

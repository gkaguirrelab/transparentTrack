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
[reconstructedPupilAzi, reconstructedPupilEle, reconstructedPupilArea] = ...
    pupilProjection_inv(pupilEllipseOnImagePlane, eyeCenterOfRotation, eyeRadius, projectionModel);

% Calculate the projected ellipse on the image plane for the reconstructed
% azimuth and elevation
projectedEllipseOnImagePlane = ...
    pupilProjection_fwd(reconstructedPupilAzi, reconstructedPupilEle, reconstructedPupilArea, eyeCenterOfRotation, eyeRadius, projectionModel);

predictedX = projectedEllipseOnImagePlane(1);
predictedY = projectedEllipseOnImagePlane(2);
closestXidx = 1;
closestYidx = 1;

% Obtain the x and y position of the projection of a pupil at this
% azimuth and elevation onto the pupil plane. We do this for the passed
% scene geometry, incorporating our uncertainty in the scene geometry to
% establish the range of X and Y points that could be plausible.
% 
% predictedX = [];
% predictedY = [];
% 
% for x=1:2
%     for y=1:2
%         for z=1:2
%             for r=1:2
%                 tmp_eyeCenterOfRotation = [sceneGeometry.eyeCenter.X_bounds(x) ...
%                     sceneGeometry.eyeCenter.Y_bounds(y) ...
%                     sceneGeometry.eyeCenter.Z_bounds(z)];
%                 tmp_eyeRadius = sceneGeometry.eyeRadius_bounds(r);
%                 tmp_projectedEllipseOnImagePlane = pupilProjection_fwd(reconstructedPupilAzi, reconstructedPupilEle, nan, tmp_eyeCenterOfRotation, tmp_eyeRadius, projectionModel);
%                 predictedX=[predictedX tmp_projectedEllipseOnImagePlane(1)];
%                 predictedY=[predictedY tmp_projectedEllipseOnImagePlane(2)];
%             end
%         end
%     end
% end
% 
% [~,closestXidx] = min(abs(predictedX - pupilEllipseOnImagePlane(1)));
% [~,closestYidx] = min(abs(predictedY - pupilEllipseOnImagePlane(2)));

% First constraint
%  Ceq reflects the Euclidean distance between the predicted and passed
%  center of the ellipse
c = constraintFactor .* sqrt( (pupilEllipseOnImagePlane(1) - predictedX(closestXidx)).^2 + (pupilEllipseOnImagePlane(2) - predictedY(closestYidx)).^2  );

% Second constraint
%  unused
ceq = [];

end

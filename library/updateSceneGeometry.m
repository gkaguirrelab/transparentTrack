function sceneGeometryOut = updateSceneGeometry( sceneGeometryIn, x )
% Update components of the sceneGeometry that impact gaze and glint
%
% Syntax:
%   sceneGeometryOut = updateSceneGeometry( sceneGeometryIn, x )
%
% Description:
%   Components of the sceneGeometry define the appearance of the pupil and
%   glint. This function receives a vector of parameters that are then
%   used to update properties of the sceneGeometry. This is useful in the
%   setting of a search across these parameters.
%
% Inputs:
%   sceneGeometryIn       - Structure. See createSceneGeometry.m
%   x                     - 1x8 vector. The elements x are:
%                            x(1) - camera torsion
%                            x(2:4) - camera position
%                            x(5) - joint eye rotation center scaler
%                            x(6) - differential eye rotation center scaler
%                            x(7) - joint corneal curvature scaler
%                            x(8) - differential corneal curvature scaler
%
% Outputs:
%   sceneGeometryOut     - Structure. See createSceneGeometry.m
%


% Copy the sceneGeometry from input to output
sceneGeometryOut = sceneGeometryIn;

% Store the camera torsion
sceneGeometryOut.cameraPosition.torsion = x(1);

% Store the extrinsic camera translation vector
sceneGeometryOut.cameraPosition.translation = x(2:4)';

% Extract the eye field
eye = sceneGeometryIn.eye;

% Obtain the default rotation centers for this eye
rotationCenters = human.rotationCenters( eye );

% Scale the rotation center values by the joint and differential
% parameters
sceneGeometryOut.eye.rotationCenters.azi = rotationCenters.azi .* x(5) .* x(6);
sceneGeometryOut.eye.rotationCenters.ele = rotationCenters.ele .* x(5) ./ x(6);

% Create the updated cornea 
radii = quadric.radii(eye.cornea.front.S);
D = @(radius) (radii(1) * 337.5) ./ radius.^2;
kvals = D(radii(2:3)) .* x(7);
kvals(1) = kvals(1) * x(8);
kvals(2) = kvals(2) / x(8);
eye.meta.measuredCornealCurvature = kvals;
cornea = human.cornea( eye );

% Update the glint optical system with the new tearfilm
tearFilmIdx = strcmp(sceneGeometryIn.refraction.glint.surfaceLabels,'cornea.tearfilm');
opticalSystemGlint = sceneGeometryIn.refraction.glint.opticalSystem;
opticalSystemGlint(tearFilmIdx,1:10) = cornea.tears.S;
sceneGeometryOut.refraction.glint.opticalSystem = opticalSystemGlint;

% Update the stopToCamera optical system with the new cornea
corneaBackIdx = find(strcmp(sceneGeometryIn.refraction.stopToCamera.surfaceLabels,'cornea.back'));
opticalSystemStopToCamera = sceneGeometryIn.refraction.stopToCamera.opticalSystem;
opticalSystemStopToCamera(corneaBackIdx:corneaBackIdx+2,1:10) = cornea.S;
sceneGeometryOut.refraction.stopToCamera.opticalSystem = opticalSystemStopToCamera;


end

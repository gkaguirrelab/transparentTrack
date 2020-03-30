function updateScene( obj, x )
% Update components of the sceneGeometry that impact gaze and glint
%
% Syntax:
%   obj.updateScene( x )
%
% Description:
%   Components of the sceneGeometry define the appearance of the pupil and
%   glint. This function receives a vector of parameters that are then
%   used to update properties of the sceneGeometry. This is useful in the
%   setting of a search across these parameters.
%
% Inputs:
%   sceneGeometryIn       - Structure. See createSceneGeometry.m
%   x                     - 1x13 vector. The elements x are:
%                            x(1:5) - kvals
%                            x(6:7) - [joint, differential] eye rotation scaler
%                            x(8:9) - primary eye position [azi ele]
%                            x(10) - camera torsion
%                            x(11:13) - camera translation
%
% Outputs:
%   sceneGeometryOut     - Structure. See createSceneGeometry.m
%


% Retrieve the current sceneGeometry from the object
sceneGeometryIn = obj.sceneGeometry;

% Copy the sceneGeometry from input to output
sceneGeometryOut = sceneGeometryIn;

% Extract the eye field
eye = sceneGeometryIn.eye;

% Update the eye meta data with the values from x
eye.meta.kvals = x(1:5);
eye.meta.rotationCenterScalers = x(6:7);
eye.meta.primaryPosition = x(8:9);

% Update the rotation centers
rotationCenters = human.rotationCenters( eye );
sceneGeometryOut.eye.rotationCenters = rotationCenters;

% Update the cornea
cornea = human.cornea( eye );
sceneGeometryOut.eye.cornea = cornea;

% Update the glint optical system with the new tearfilm
tearFilmIdx = strcmp(sceneGeometryIn.refraction.glint.surfaceLabels,'cornea.tearfilm');
opticalSystemGlint = sceneGeometryIn.refraction.glint.opticalSystem;
opticalSystemGlint(tearFilmIdx,1:10) = cornea.tears.S;
sceneGeometryOut.refraction.glint.opticalSystem = opticalSystemGlint;

% Update the stopToMedium optical system with the new cornea
corneaBackIdx = find(strcmp(sceneGeometryIn.refraction.stopToMedium.surfaceLabels,'cornea.back'));
opticalSystemStopToMedium = sceneGeometryIn.refraction.stopToMedium.opticalSystem;
opticalSystemStopToMedium(corneaBackIdx:corneaBackIdx+2,1:10) = cornea.S;
sceneGeometryOut.refraction.stopToMedium.opticalSystem = opticalSystemStopToMedium;

% Store the camera torsion
sceneGeometryOut.cameraPosition.torsion = x(10);

% Store the extrinsic camera translation vector
sceneGeometryOut.cameraPosition.translation = x(11:13)';

% Store the updated sceneGeometry in the object
obj.sceneGeometry = sceneGeometryOut;

end

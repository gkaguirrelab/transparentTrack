function updateScene( obj )
% Update components of the sceneGeometry that impact gaze and glint
%
% Syntax:
%  obj.updateScene()
%
% Description:
%   Components of the sceneGeometry define the appearance of the pupil and
%   glint. This function receives a vector of parameters that are then
%   used to update properties of the sceneGeometry. This is useful in the
%   setting of a search across these parameters.
%
% Inputs:
%   none. All contained in the object.
%
% Outputs:
%   sceneGeometryOut     - Structure. See createSceneGeometry.m
%


% Retrieve the current sceneGeometry from the object
sceneGeometryIn = obj.sceneGeometry;

% Copy the sceneGeometry from input to output
sceneGeometryOut = sceneGeometryIn;

% Get the model parameters and model description
x = obj.x;
model = obj.model;

% Extract the eye field
eye = sceneGeometryIn.eye;

% Update the eye meta data with the values from x
eye.meta.kvals = x(model.func.fieldSetIdx('eye','kvals'));
eye.meta.rotationCenterScalers = x(model.func.fieldSetIdx('eye','rotationCenterScalers'));
eye.meta.primaryPosition = x(model.func.fieldSetIdx('scene','primaryPosition'));

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
sceneGeometryOut.cameraPosition.torsion = x(model.func.fieldParamIdx('scene','torsion'));

% Store the extrinsic camera translation vector
sceneGeometryOut.cameraPosition.translation = x(model.func.fieldSetIdx('scene','translation'))';

% Store the updated sceneGeometry in the object
obj.sceneGeometry = sceneGeometryOut;

end

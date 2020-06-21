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

% Get the setup args
setupArgs = obj.setupArgs;

% Retrieve the current sceneGeometry from the object
sceneGeometryIn = obj.sceneGeometry;

% Copy the sceneGeometry from input to output
sceneGeometryOut = sceneGeometryIn;

% Get the model parameters and model description
x = obj.x;
xLast = obj.xLast;
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

% Only re-calculate the optical systems if there has been a change in the
% kvals
if any( x(model.func.fieldSetIdx('eye','kvals')) ~= xLast(model.func.fieldSetIdx('eye','kvals')) )
    
    % Update the glint optical system
    sceneGeometryOut.refraction.glint.opticalSystem = ...
        assembleOpticalSystem( sceneGeometryOut.eye, 'surfaceSetName', 'glint', 'skipMagCalc', true, setupArgs{:});
    
    % Update the stopToMedium optical system with the cornea
    sceneGeometryOut.refraction.stopToMedium.opticalSystem = ...
        assembleOpticalSystem( sceneGeometryOut.eye, 'surfaceSetName', 'stopToMedium', 'skipMagCalc', true, setupArgs{:});
    
end

% Store the camera torsion
sceneGeometryOut.cameraPosition.torsion = x(model.func.fieldParamIdx('scene','torsion'));

% Store the extrinsic camera translation vector
sceneGeometryOut.cameraPosition.translation = x(model.func.fieldSetIdx('scene','translation'))';

% Add any common depth adjustments acros scenes
sceneGeometryOut.cameraPosition.translation(3) = ...
    sceneGeometryOut.cameraPosition.translation(3) + ...
    x(model.func.fieldSetIdx('eye','commonDepth'));

% Store the updated sceneGeometry in the object
obj.sceneGeometry = sceneGeometryOut;

end

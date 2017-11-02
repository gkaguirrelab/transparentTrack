function [eccentricity, theta] = constrainEllipseBySceneGeometry (ellipseCenter,sceneGeometry, varargin)
% [eccentricity, theta] = constrainEllipseBySceneGeometry (ellipseCenter,sceneGeometry)
%
% This function returns the expected ecceentricity and tilt for an ellipse,
% given the location of the center and the scene geometry.
% The function can either use the Z coordinate for the eyeball in the scene
% geometry as the distance of the eye from the scene, or a user input
% distance range in pixel (2 elements vector). In the latter case, the
% routine will return the range of expected eccentricity and a single theta.
% 
% Output
% eccentricity - either a one element or 2 element vector with the expected
%   eccentricity value or range.
% theta -the expected tilt value for the ellipse (does not depend on the distance).
% 
% Input (required)
% ellipseCenter - [X Y] coordinate for the ellipse center. Note that
%   passing the full parametrization of a transparent ellipse will work as
%   well.
% sceneGeometry - struct with scene geometry
% 
% Optional Input (analysis)
% distanceFromSceneRangePx - [distInPx] or [minDistPx maxDistPx] estimate
%   of the distance of the eye from the scene plane (i.e. from the camera).
%   This must be in pixel units.


%% input parser
p = inputParser; p.KeepUnmatched = true;

% required input
p.addRequired('ellipseCenter',@isnumeric);
p.addRequired('sceneGeometry',@isstruct);

% optional input
p.addParameter('distanceFromScenePx', [],@isnumeric)

%parse
p.parse(ellipseCenter,sceneGeometry, varargin{:})


%% derive rotation angles and reconstruct eccentricity and theta

    % derive pupil center cartesian 3D coordinates
    pupilCenter.X = ellipseCenter(1);
    pupilCenter.Y = ellipseCenter(2);
    pupilCenter.Z = sqrt(sceneGeometry.eyeballRadius^2 - (pupilCenter.X - sceneGeometry.eyeballCenter.X)^2 - (pupilCenter.Y - sceneGeometry.eyeballCenter.Y)^2) +sceneGeometry.eyeballCenter.Z ;
    pupilAzi = atand((pupilCenter.X - sceneGeometry.eyeballCenter.X)/(pupilCenter.Z - sceneGeometry.eyeballCenter.Z));
    pupilEle = atand((pupilCenter.Y - sceneGeometry.eyeballCenter.Y)/(pupilCenter.Z - sceneGeometry.eyeballCenter.Z));
    
    reconstructedTransparentEllipse = pupilProjection_fwd(pupilAzi, pupilEle, [sceneGeometry.eyeballCenter.X sceneGeometry.eyeballCenter.Y sceneGeometry.eyeballCenter.Z]);
    eccentricity = reconstructedTransparentEllipse(4);
    theta = reconstructedTransparentEllipse(5);
    
end

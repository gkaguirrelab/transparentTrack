function [eccentricity, theta] = constrainEllipseBySceneGeometry (ellipseCenter, sceneGeometry, varargin)
% [eccentricity, theta] = constrainEllipseBySceneGeometry (ellipseCenter,sceneGeometry)
%
% This function returns the expected ecceentricity and tilt for an ellipse,
% given the location of the center and the scene geometry.
% The function can either use the Z coordinate for the eye in the scene
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


%% input parser
p = inputParser; p.KeepUnmatched = true;

% required input
p.addRequired('ellipseCenter',@isnumeric);
p.addRequired('sceneGeometry',@isstruct);


%parse
p.parse(ellipseCenter,sceneGeometry, varargin{:})


%% derive rotation angles and reconstruct eccentricity and theta

    % derive pupil center cartesian 3D coordinates
    pupilCenter.X = ellipseCenter(1);
    pupilCenter.Y = ellipseCenter(2);
    pupilCenter.Z = sqrt(sceneGeometry.eyeRadius^2 - (pupilCenter.X - sceneGeometry.eyeCenter.X)^2 - (pupilCenter.Y - sceneGeometry.eyeCenter.Y)^2) +sceneGeometry.eyeCenter.Z ;
    pupilAzi = atand((pupilCenter.X - sceneGeometry.eyeCenter.X)/(pupilCenter.Z - sceneGeometry.eyeCenter.Z));
    pupilEle = atand((pupilCenter.Y - sceneGeometry.eyeCenter.Y)/(pupilCenter.Z - sceneGeometry.eyeCenter.Z));
    
    reconstructedTransparentEllipse = pupilProjection_fwd(pupilAzi, pupilEle, [sceneGeometry.eyeCenter.X sceneGeometry.eyeCenter.Y sceneGeometry.eyeCenter.Z]);
    eccentricity = reconstructedTransparentEllipse(4);
    theta = reconstructedTransparentEllipse(5);
    
end

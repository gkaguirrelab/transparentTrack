function [eccentricity, theta] = constrainEllipseBySceneGeometry (ellipseCenter,sceneGeometryFileName, varargin)
% [eccentricity, theta] = constrainEllipseBySceneGeometry (ellipseCenter,sceneGeometryFileName)
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
% sceneGeometryFileName - name of the scene geometry file
% 
% Optional Input (analysis)
% distanceFromSceneRangePx - [minDistPx maxDistPx] estimate of the distance
%   of the eye from the scene plane (i.e. from the camera). This must be in
%   pixel units.


%% input parser
p = inputParser; p.KeepUnmatched = true;

% required input
p.addRequired('ellipseCenter',@isnumeric);
p.addRequired('sceneGeometryFileName',@isstr);

% optional input
p.addParameter('distanceFromSceneRangePx', [],@isnumeric)

%parse
p.parse(ellipseCenter, sceneGeometryFileName, varargin{:})
%% load scene geometry
load(sceneGeometryFileName)

%% derive rotation angles and reconstruct eccentricity and theta
if isempty (p.Results.distanceFromSceneRangePx)
    pupilAzi = atand((ellipseCenter(1)-sceneGeometry.eyeball.X)/sceneGeometry.eyeball.Z);
    pupilEle = atand((ellipseCenter(2)-sceneGeometry.eyeball.Y)/sceneGeometry.eyeball.Z);
    
    reconstructedTransparentEllipse = pupilProjection_fwd(pupilAzi, pupilEle, [sceneGeometry.eyeball.X sceneGeometry.eyeball.Y sceneGeometry.eyeball.Z]);
    eccentricity = reconstructedTransparentEllipse(4);
    theta = reconstructedTransparentEllipse(5);
else
    % calculate lower bound values for eccentricity
    minPupilAzi = atand((ellipseCenter(1)-sceneGeometry.eyeball.X)/p.Results.distanceFromSceneRangePx(1));
    minPupilEle = atand((ellipseCenter(2)-sceneGeometry.eyeball.Y)/p.Results.distanceFromSceneRangePx(1));
    minReconstructedTransparentEllipse = pupilProjection_fwd(minPupilAzi, minPupilEle, [sceneGeometry.eyeball.X sceneGeometry.eyeball.Y p.Results.distanceFromSceneRangePx(1)]);
    eccentricity(1) = minReconstructedTransparentEllipse(4);
    
    % calculate upper bound values fpr eccentricity
    maxPupilAzi = atand((ellipseCenter(1)-sceneGeometry.eyeball.X)/p.Results.distanceFromSceneRangePx(2));
    maxPupilEle = atand((ellipseCenter(2)-sceneGeometry.eyeball.Y)/p.Results.distanceFromSceneRangePx(2));
    maxReconstructedTransparentEllipse = pupilProjection_fwd(maxPupilAzi, maxPupilEle, [sceneGeometry.eyeball.X sceneGeometry.eyeball.Y p.Results.distanceFromSceneRangePx(2)]);
    eccentricity(2) = maxReconstructedTransparentEllipse(4);
    
    % calculate the tilt
    theta = maxReconstructedTransparentEllipse(5);
end

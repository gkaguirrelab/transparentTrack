function [reconstructedPupilAzi, reconstructedPupilEle, reconstructedPupilArea] = pupilProjection_inv(transparentEllipse, eyeCenter, eyeRadius, projectionModel)
% [reconstructedHorizontalAngle, reconstructedVerticalAngle] =
% pupilProjection_inv(transparentEllipse,centerOfProjection)
%
% Returns the horiziontal and vertical angles of tilt of the pupil plane
% with reference to the scene plane using the params of the transparent
% ellipse and the coordinates of the center of projection on the scene.
%
% Note that the linear units must be uniform (eg. all pixels or all mm) for
% both the transparent ellipse (where applicable) and the center of
% projection coordinates.
%
% ref http://web.ncf.ca/aa456/scale/ellipse.html
%
% Outputs:
%   reconstructedPupilAzi - rotation of the pupil in the XY plane in
%       degrees, with the center being the centerOfProjection on the scene
%   reconstructedPupilEle - elevation of the pupil from the XY plane in
%       degrees, with the center being the centerOfProjection on the scene.
%   reconstructedPupilArea - pupil area, in the same (squared) units as the
%       transparent ellipse center.
%
% Required inputs:
%   transparentEllipse - ellipse in transparent form
%   eyeCenter - 3D coordinates of the eye center in the scene reference
%       system. This could be a vector assembled from the .X, .Y, and .Z
%       vaues of the sceneGeometry.eyeCenter
%   eyeRadius - the estimate radius for the eye (can be empty for
%       orthogonal projection)
%   projectionModel - string that identifies the projection model to use.
%       Options include "orthogonal" and "pseudoPerspective"


%% variable definitions

% if we have a non-defined ellipse, just return all nans
if any(isnan(transparentEllipse(1:2)))
    reconstructedPupilAzi = nan;
    reconstructedPupilEle = nan;
    reconstructedPupilArea = nan;
    return
end

% projection variables
k = sqrt(1 - (transparentEllipse(4)^2));
theta = (transparentEllipse(5));

% keep theta between zero and pi
if theta < 0
    theta = theta + pi;
end

% pupil center
switch projectionModel
    case 'orthogonal'
        centerX = (transparentEllipse(1));
        centerY = (transparentEllipse(2));
    case 'pseudoPerspective'
        % at this stage we just use the orthogonal projection of the
        % ellipse center to determine the location of the ellipse with
        % respect to the center of projection. we will correct this below.
        centerX = (transparentEllipse(1));
        centerY = (transparentEllipse(2));
    otherwise
        error('I do not know that projection case');
end

% derive horizontal tilt angle (azimuth)
if ~any(isnan(eyeCenter))
    % If we have an eyeCenter defined, then we can provide a unique
    % solution for the azimuth and elevation
    if centerX > eyeCenter(1)
        reconstructedPupilAzi = asind(sqrt((sin(theta))^2 * (1 -k^2)));
    elseif centerX < eyeCenter(1)
        reconstructedPupilAzi = -asind(sqrt((sin(theta))^2 * (1 -k^2)));
    elseif centerX == eyeCenter(1)
        reconstructedPupilAzi = 0;
    end
    
    % derive vertical tilt angle (elevation)
    if centerY > eyeCenter(2)
        reconstructedPupilEle = asind(sqrt(((cos(theta))^2 *(k^2 - 1))/(((sin(theta))^2 * (1 -k^2)) -1)));
    elseif centerY < eyeCenter(2)
        reconstructedPupilEle = -asind(sqrt(((cos(theta))^2 *(k^2 - 1))/(((sin(theta))^2 * (1 -k^2)) -1)));
    elseif centerY == eyeCenter(2)
        reconstructedPupilEle = 0;
    end
else
    % If we do not have a defined eyeCenter, then we return both of the
    % possible solutions
    if theta < pi/2
        reconstructedPupilAzi(1) = asind(sqrt((sin(theta))^2 * (1 -k^2)));
        reconstructedPupilAzi(2) = -asind(sqrt((sin(theta))^2 * (1 -k^2)));
        reconstructedPupilEle(1) = -asind(sqrt(((cos(theta))^2 *(k^2 - 1))/(((sin(theta))^2 * (1 -k^2)) -1)));
        reconstructedPupilEle(2) = asind(sqrt(((cos(theta))^2 *(k^2 - 1))/(((sin(theta))^2 * (1 -k^2)) -1)));
    else
        reconstructedPupilAzi(1) = asind(sqrt((sin(theta))^2 * (1 -k^2)));
        reconstructedPupilAzi(2) = -asind(sqrt((sin(theta))^2 * (1 -k^2)));
        reconstructedPupilEle(1) = asind(sqrt(((cos(theta))^2 *(k^2 - 1))/(((sin(theta))^2 * (1 -k^2)) -1)));
        reconstructedPupilEle(2) = -asind(sqrt(((cos(theta))^2 *(k^2 - 1))/(((sin(theta))^2 * (1 -k^2)) -1)));
    end
end


% If we are given an ellipse area, calculate a pupil area
if ~isnan(transparentEllipse(3))
    switch projectionModel
        case 'orthogonal'
            reconstructedPupilRadius = sqrt(transparentEllipse(3) / (pi * k));
            reconstructedPupilArea = pi * reconstructedPupilRadius^2;
        case 'pseudoPerspective'
            if ~any(isnan(eyeCenter))
                
                % calculate the perspective correction factor
                sceneDistance = abs(eyeCenter(3) - eyeRadius);
                relativeDepth = eyeRadius*(1-(cosd(reconstructedPupilEle)*cosd(reconstructedPupilAzi)));
                perspectiveCorrectionFactor = (sceneDistance/(sceneDistance + relativeDepth));
                
                % calculate pupil radius including the perspective
                % correction factor for the ellipse area
                apparentPupilRadius = sqrt(transparentEllipse(3)/ (pi * k));                
                reconstructedPupilRadius = apparentPupilRadius./perspectiveCorrectionFactor;
                
                % calculate the area
                reconstructedPupilArea = pi * reconstructedPupilRadius^2;
            else
                error('Cannot apply perspective correction without knowing the center of projection')
            end
    end
else
    reconstructedPupilArea = nan;
end

end % function

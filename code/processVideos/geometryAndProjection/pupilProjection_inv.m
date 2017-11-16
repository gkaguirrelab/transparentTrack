function [reconstructedPupilAzi, reconstructedPupilEle, reconstructedPupilArea] = pupilProjection_inv(transparentEllipse, centerOfProjection, eyeRadius, sceneDistance, projectionModel)
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
%   center of projection -[X Y] coordinates of the center of projection on
%       the scene. The center of projection is the center of the projected
%       pupil that preserves its circular shape in the scene.
%   eyeRadius - the estimate radius for the eye
%   sceneDistance - orthogonal distance of the scene from the center of
%       projection
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
theta = round(transparentEllipse(5),4);

% pupil center
switch projectionModel
    case 'orthogonal'
        centerX = round(transparentEllipse(1),4);
        centerY = round(transparentEllipse(2),4);
    case 'pseudoPerspective'
        % at this stage we just use the ortogonal projection of the ellipse
        % center to determine the location of the ellipse with respect to
        % the center of correction.
        centerX = round(transparentEllipse(1),4);
        centerY = round(transparentEllipse(2),4);
end

% derive horizontal tilt angle (azimuth)
if ~isnan(centerOfProjection)
    if centerX > centerOfProjection(1)
        reconstructedPupilAzi = asind(sqrt((sin(theta))^2 * (1 -k^2)));
    elseif centerX < centerOfProjection(1)
        reconstructedPupilAzi = - asind(sqrt((sin(theta))^2 * (1 -k^2)));
    elseif centerX == centerOfProjection(1)
        reconstructedPupilAzi = 0;
    end
    
    % derive vertical tilt angle (elevation)
    if centerY > centerOfProjection(2)
        reconstructedPupilEle = asind(sqrt(((cos(theta))^2 *(k^2 - 1))/(((sin(theta))^2 * (1 -k^2)) -1)));
    elseif centerY < centerOfProjection(2)
        reconstructedPupilEle = -asind(sqrt(((cos(theta))^2 *(k^2 - 1))/(((sin(theta))^2 * (1 -k^2)) -1)));
    elseif centerY == centerOfProjection(2)
        reconstructedPupilEle = 0;
    end
else
    reconstructedPupilAzi(1) = asind(sqrt((sin(theta))^2 * (1 -k^2)));
    reconstructedPupilAzi(2) = -asind(sqrt((sin(theta))^2 * (1 -k^2)));
    reconstructedPupilEle(1) = asind(sqrt(((cos(theta))^2 *(k^2 - 1))/(((sin(theta))^2 * (1 -k^2)) -1)));
    reconstructedPupilEle(2) = - asind(sqrt(((cos(theta))^2 *(k^2 - 1))/(((sin(theta))^2 * (1 -k^2)) -1)));
end


% see if you can derive the pupil radius
if ~isnan(transparentEllipse(3))
    switch projectionModel
        case 'orthogonal'
            reconstructedPupilRadius = sqrt(transparentEllipse(3) / (pi * k));
            reconstructedPupilArea = pi * reconstructedPupilRadius^2;
        case 'pseudoPerspective'
            if ~isnan(centerOfProjection)
                % calculate the perspective correction factor
                relativeDepth = eyeRadius*(1-(cosd(reconstructedPupilEle)*cosd(reconstructedPupilAzi)));
                perspectiveCorrectionFactor = (sceneDistance/(sceneDistance + relativeDepth));
                
                % calculate pupil radius including the perspective
                % correction factor for the ellipse area
                reconstructedPupilRadius = sqrt(transparentEllipse(3)*perspectiveCorrectionFactor / (pi * k));
                
                % calculate the area
                reconstructedPupilArea = pi * reconstructedPupilRadius^2;    
            else
                error('Cannot apply perspective correction without knowing the center of projection')
            end
    end
end

end % function

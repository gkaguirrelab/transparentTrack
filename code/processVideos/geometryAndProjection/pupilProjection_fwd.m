function reconstructedTransparentEllipse = pupilProjection_fwd(pupilAzi, pupilEle, pupilArea, eyeCenter, eyeRadius, projectionModel)
% reconstructedTransparentEllipse = pupilProjection_fwd(pupilAzi, pupilEle, pupilCenter3D)
%
% Returns the transparent ellipse params of the pupil projection on the
% scene, using the 3D pupil center coordinates and the horiziontal and
% vertical angles of tilt (in degrees) of the pupil center with reference
% to the scene plane.
%
% If the pupil area is available, the transparent param for the ellipse
% area will return the projected ellipse area, otherwise it will be left as
% NaN.
%
% Note that the linear units must be uniform (eg. all pixels or all mm) for
% both the pupil center and the transparent ellipse parameters (where
% applicable).
%
% Outputs:
%   reconstructedTransparentEllipse - ellipse in transparent form
%
% Required inputs:
%   pupilAzi - rotation of the pupil in the XY plane in degrees, with
%       the center being the centerOfProjection on the scene.
%   pupilEle - elevation of the pupil from the XY plane in degrees, with
%       the center being the centerOfProjection on the scene.
%   pupilArea - if not set to nan, the routine will return the ellipse area
%   eyeCenter - 3D coordinates of the eye center in the scene reference
%       system. This could be a vector assembled from the .X, .Y, and .Z
%       vaues of the sceneGeometry.eyeCenter
%   eyeRadius - the estimate radius for the eye (can be empty for
%       orthogonal projection)
%   projectionModel - string that identifies the projection model to use.
%       Options include "orthogonal" and "pseudoPerspective"


%% main

% initiate reconstructedTransparentEllipse
reconstructedTransparentEllipse = nan(1,5);

% if we have a non-defined case, just return all nans
if any(isnan([pupilAzi pupilEle]))
    return
end

% calculate the pupilCenter3D
pupilCenter3D(1) = eyeRadius*(cosd(pupilEle)*sind(pupilAzi));
pupilCenter3D(2) = eyeRadius*sind(pupilEle);
pupilCenter3D(3) = eyeRadius*(cosd(pupilEle)*cosd(pupilAzi));

% define ellipse center
switch projectionModel
    case 'orthogonal'
        % under orthogonal hypothesis the ellipse center in 2D is coincident with
        % the ellipse center in the plane of projection.
        reconstructedTransparentEllipse(1) = pupilCenter3D(1);
        reconstructedTransparentEllipse(2) = pupilCenter3D(2);
    case 'pseudoPerspective'
        % for the pseudoPerspective correction, we uniformly scale the
        % orthogonal projection according to the scene distance and the eye
        % radius.
        
        % get the perspective projection correction factor
        relativeDepth = eyeRadius*(1-(cosd(pupilEle)*cosd(pupilAzi)));
        sceneDistance = abs(eyeCenter(3) - eyeRadius);
        perspectiveCorrectionFactor = sceneDistance/(sceneDistance + relativeDepth);
        
        % apply perspective correction factor to the pupil center
        reconstructedTransparentEllipse(1) = pupilCenter3D(1) * perspectiveCorrectionFactor;
        reconstructedTransparentEllipse(2) = pupilCenter3D(2) * perspectiveCorrectionFactor;
end

% derive transparent parameters
% eccentricity
e = sqrt((sind(pupilEle))^2 - ((sind(pupilAzi))^2 * ((sind(pupilEle))^2 - 1)));
reconstructedTransparentEllipse(4) = e;

theta = nan;
% tilt
if pupilAzi > 0  &&  pupilEle > 0
    theta = - asin(sind(pupilAzi)/e);
elseif pupilAzi > 0 && pupilEle < 0
    theta =  asin(sind(pupilAzi)/e);
elseif pupilAzi < 0  &&  pupilEle > 0
    theta = pi/2 + asin(sind(pupilAzi)/e);
elseif pupilAzi < 0 && pupilEle < 0
    theta =  pi/2 - asin(sind(pupilAzi)/e);
elseif pupilAzi == 0 && pupilEle == 0
    theta = 0;
elseif pupilAzi ==  0 && pupilEle ~= 0
    theta = 0;
elseif pupilEle == 0 && pupilAzi ~= 0
    theta = pi/2;
else
    % Couldn't constrain the theta
    % TO INVESTIGATE: If it is the case that this occurs when the ellipse
    % if very close to circular, we may choose to set theta to zero in this
    % event.
    theta = nan;
end
reconstructedTransparentEllipse(5) = theta;

% area (if pupilRadius available)
if ~isnan(pupilArea)
    switch projectionModel
        case 'orthogonal'
            % under orthogonal hypothesis, the semimajor axis of the
            % ellipse equals the pupil circle radius. If that is known, it
            % can be assigned and used later to determine the ellipse area.
            pupilRadius = sqrt(pupilArea/pi);
            semiMajorAxis = pupilRadius;
        case 'pseudoPerspective'
            % apply scaling factor to pupil area
            pupilRadius = sqrt(pupilArea/pi);
            semiMajorAxis = pupilRadius.*perspectiveCorrectionFactor;
    end
    semiMinorAxis = semiMajorAxis*sqrt(1-e^2);
    reconstructedTransparentEllipse(3) = pi*semiMajorAxis*semiMinorAxis;
end

end % function


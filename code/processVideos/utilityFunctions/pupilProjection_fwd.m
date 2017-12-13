function reconstructedTransparentEllipse = pupilProjection_fwd(pupilAzi, pupilEle, pupilArea, eyeCenter, eyeRadius, projectionModel)
% Transform pupil azimuth, elevation, area to an ellipse on the image plane
%
% Description:
%	Returns the transparent ellipse params of the pupil projection on the
%   image plane, using the horiziontal and vertical angles of tilt (in
%   degrees) of the pupil center along with the scene geometry.
%
%   Note that we use degrees to specify pupil azimuth and elevation, but
%   radians for the theta value in the transparent ellipse formulation.
%   This is in part to help us keep the two units separate conceptually.
%
% Inputs:
%   pupilAzi    - Rotation of the pupil in the XY plane in degrees, with 
%                 the center being the centerOfProjection on the scene
%   pupilEle    - Elevation of the pupil from the XY plane in degrees, with
%                 the center being the centerOfProjection on the scene              
%   pupilArea   - Pupil area, in the same (squared) units as the 
%                 transparent ellipse center
%   eyeCenter   - 3D coordinates of the eye center in the scene reference
%                 system. This could be a vector assembled from the .X, .Y,
%                 and .Z vaues of the sceneGeometry.eyeCenter
%   eyeRadius   - the radius of the eye (this variable may be empty for
%                 the orthogonal projection case)
%   projectionModel - string that identifies the projection model to use.
%                 Options include "orthogonal" and "pseudoPerspective"
%
% Outputs:
%   reconstructedTransparentEllipse  - Ellipse in transparent form
%


% initiate the reconstructedTransparentEllipse
reconstructedTransparentEllipse = nan(1,5);

% if we have a non-defined case, just return all nans
if any(isnan([pupilAzi pupilEle]))
    return
end

% calculate the pupilCenter3D relative to the eyeCenter
pupilCenter3D(1) = eyeRadius*sind(pupilAzi)*cosd(pupilEle);
pupilCenter3D(2) = eyeRadius*sind(pupilEle);
pupilCenter3D(3) = eyeRadius - eyeRadius*cosd(pupilAzi)*cosd(pupilEle);

% define ellipse center
switch projectionModel
    case 'orthogonal'
        % under orthogonal hypothesis the ellipse center in 2D is
        % coincident with the ellipse center in the plane of projection.
        reconstructedTransparentEllipse(1) = pupilCenter3D(1) + eyeCenter(1);
        reconstructedTransparentEllipse(2) = pupilCenter3D(2) + eyeCenter(2);
    case 'pseudoPerspective'
        % for the pseudoPerspective correction, we uniformly scale the
        % orthogonal projection according to the scene distance and the eye
        % radius.
        
        % get the perspective projection correction factor
        sceneDistance = eyeCenter(3) - eyeRadius;
        perspectiveCorrectionFactor = sceneDistance/(sceneDistance + pupilCenter3D(3));
        
        % apply perspective correction factor to the pupil center
        reconstructedTransparentEllipse(1) = (pupilCenter3D(1) * perspectiveCorrectionFactor) + eyeCenter(1);
        reconstructedTransparentEllipse(2) = (pupilCenter3D(2) * perspectiveCorrectionFactor) + eyeCenter(2);
end

% derive ellipse eccentricity
e = sqrt((sind(pupilEle))^2 - ((sind(pupilAzi))^2 * ((sind(pupilEle))^2 - 1)));
reconstructedTransparentEllipse(4) = e;

% derive ellipse theta
theta = nan;
if pupilAzi > 0  &&  pupilEle > 0
    theta = - asin(sind(pupilAzi)/e);
elseif pupilAzi > 0 && pupilEle < 0
    theta =  asin(sind(pupilAzi)/e);
elseif pupilAzi < 0  &&  pupilEle > 0
    theta = - asin(sind(pupilAzi)/e);
elseif pupilAzi < 0 && pupilEle < 0
    theta =  asin(sind(pupilAzi)/e);
elseif abs(pupilAzi) < 1e-12 && abs(pupilEle) < 1e-12
    theta = 0;
elseif abs(pupilAzi) < 1e-12 && abs(pupilEle) >= 1e-12
    theta = 0;
elseif abs(pupilEle) < 1e-12 && abs(pupilAzi) >= 1e-12
    theta = pi/2;
else
    % Couldn't constrain the theta. This shouldn't happen
    warning('For some reason the theta was unconstrained. Leaving as nan'); 
end

% Keep the theta values between 0 and pi
if theta < 0
    theta = theta+pi;
end
if theta > pi
    theta = theta - pi;
end

reconstructedTransparentEllipse(5) = theta;

% area (if pupilArea available)
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


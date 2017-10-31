function reconstructedTransparentEllipse = pupilProjection_fwd(pupilAzi, pupilEle, pupilCenter3D, varargin)
% reconstructedTransparentEllipse = pupilProjection_fwd(pupilAzi, pupilEle, pupilCenter3D)
% 
% Returns the transparent ellipse params of the pupil projection on the
% scene, using the 3D pupil center coordinates and the horiziontal and
% vertical angles of tilt (in degrees) of the pupil center with reference
% to the scene plane.
% If the pupil radius is available, the transparent param for the ellipse
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
% pupilAzi - rotation of the pupil in the XY plane in degrees,
%       with the center being the centerOfProjection on the scene.
% pupilEle - elevation of the pupil from the XY plane in degrees,
%       with the center being the centerOfProjection on the scene.
% pupilCenter3D - center of the pupil in 3D, where X and Y are parallel to
%   the scene plane and Z is the distance of the center of the pupil from the
%   scene plane.
% 
% Optional inputs:
%   pupilRadius -  when available, the routine will return the ellipse area
%   perspectiveCorrection - DEV PLACEHOLDER param to appy perspective
%       correction

%% parse input and define variables
p = inputParser; p.KeepUnmatched = true;

% required input
p.addRequired('pupilAzi',@isnumeric);
p.addRequired('pupilEle',@isnumeric);
p.addRequired('pupilCenter3D',@isnumeric);
p.addRequired('centerOfProjection',@isnumeric);

% optional analysis params
p.addParameter('pupilRadius', nan, @isnumeric);
p.addParameter('perspectiveCorrection', 0, @isnumeric);

% Environment parameters
p.addParameter('tbSnapshot',[],@(x)(isempty(x) | isstruct(x)));
p.addParameter('timestamp',char(datetime('now')),@ischar);
p.addParameter('username',char(java.lang.System.getProperty('user.name')),@ischar);
p.addParameter('hostname',char(java.net.InetAddress.getLocalHost.getHostName),@ischar);

% parse
p.parse(pupilAzi, pupilEle, pupilCenter3D,varargin{:})


%% initiate reconstructedTransparentEllipse
reconstructedTransparentEllipse = nan(1,5);

%% define ellipse center
if p.Results.perspectiveCorrection == 0
    % under orthogonal hypothesis the ellipse center in 2D is coincident with
    % the ellipse center in the plane of projection.
    reconstructedTransparentEllipse(1) = pupilCenter3D(1);
    reconstructedTransparentEllipse(2) = pupilCenter3D(2);
else
    % DEV PLACEHOLDER: if there is a perspective correction value to apply,
    % apply it!
end

%% derive transparent parameters
% eccentricity
e = sqrt((sind(pupilEle))^2 - ((sind(pupilAzi))^2 * ((sind(pupilEle))^2 - 1)));
reconstructedTransparentEllipse(4) = e;

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
end
reconstructedTransparentEllipse(5) = theta;

% area (if pupilRadius available)
if ~isnan(p.Results.pupilRadius)
    if p.Results.perspectiveCorrection == 0
        % under orthogonal hypotesis, the semimajor axis of the ellipse equals the
        % pupil circle radius. If that is known, it can be assigned and used later
        % to determine the ellipse area.
        semiMajorAxis = p.Results.pupilRadius;
    elseif p.Results.perspectiveCorrection > 0
        % DEV PLACEHOLDER: if there is a perspective correction value to apply,
        % apply it!
    end
semiMinorAxis = semiMajorAxis*sqrt(1-e^2);
reconstructedTransparentEllipse(3) = pi*semiMajorAxis*semiMinorAxis;
end


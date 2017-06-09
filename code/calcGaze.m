function gaze = calcGaze(params)

% Calculate the gaze location using eye tracking data
%
%   Usage:
%       gaze = calcGaze(params)
%
%   Inputs:
%       params.pupil.X      - vector of pupil X coordinates (pixels)
%       params.pupil.Y      - vector of pupil Y coordinates (pixels)
%       params.glint.X      - vector of glint X coordinates (pixels)
%       params.glint.Y      - vector of glint Y coordinates (pixels)
%       params.viewDist     - viewing distance (mm)
%       params.rpc          - adjustment parameter (see 'calcRpc')
%       params.calMat       - transformation matrix (see 'calcCalMat');
%
%
%   Outputs:
%       gaze.X              - vector of gaze X coordinates (mm)
%       gaze.Y              - vector of gaze Y coordinates (mm)
%       gaze.ecc            - vector of gaze eccentricity values (degrees visual angle)
%       gaze.pol            - vector of gaze polar angle values (degrees)
%
%   Notes for 'gaze.pol':
%           0               - upper vertical meridian
%           90              - right horizontal meridian
%           180             - lower vertical meridian
%           270             - left horizontal meridian
%
%   Written by Andrew S Bock Oct 2016
%   Jan 2017 - updated to use liveTrackCartToPol (GF)

%% Calculate the gaze
gaze.X                  = nan(size(params.pupil.X));
gaze.Y                  = nan(size(params.pupil.X));
gaze.ecc                = nan(size(params.pupil.X));
gaze.pol                = nan(size(params.pupil.X));
tmp                     = nan(length(params.pupil.X),3);
for i = 1:length(params.pupil.X)
    pX                  = params.pupil.X(i);
    pY                  = params.pupil.Y(i);
    gX                  = params.glint.X(i);
    gY                  = params.glint.Y(i);
    aXYZW               = params.calMat * [...
        (pX-gX)/params.rpc; ...
        (pY-gY)/params.rpc; ...
        (1 - sqrt(((pX-gX)/params.rpc)^2 + ((pY-gY)/params.rpc)^2)); ...
        1];
    tmp(i,:)            = (aXYZW(1:3)/aXYZW(4))';
end
%% Pull out the X and Y values
gaze.X                  = tmp(:,1);
gaze.Y                  = tmp(:,2);
% note tmp(:,3) is the viewing distance
%% Convert to eccentricty and polar angle
for ii = 1:length(gaze.X)
    [gaze.ecc(ii), gaze.pol(ii)] = liveTrackCartToPol(gaze.X(ii),gaze.Y(ii),params.viewDist);
end
gaze.ecc = gaze.ecc';
gaze.pol = gaze.pol';


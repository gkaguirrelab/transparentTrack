function rpc = calcRpc(params)

% Function to calculate the adjustment factor needed for gaze calculation
%   
%   The glint center and pupil center are not aligned when the subject is 
%   looking straight ahead, thus an adjustment is required. The name 'rpc'
%   specifies the 'ratio of pupil to cornea'
%
%   Usage:
%       rpc = calcRpc(params)
%
%   Inputs:
%       params.pupil.X      - vector of pupil X coordinates (pixels)
%       params.pupil.Y      - vector of pupil Y coordinates (pixels)
%       params.glint.X      - vector of glint X coordinates (pixels)
%       params.glint.Y      - vector of glint Y coordinates (pixels)
%       params.targets.X    - vector of target X coordinates (pixels)
%       params.targets.Y    - vector of target Y coordinates (pixels)
%
%   Output:
%       rpc                 - scalar adjustment factor
%
%   Written by Andrew S Bock Oct 2016

%% Pull out variables
targets         = [params.targets.X, params.targets.Y];
pupil           = [params.pupil.X, params.pupil.Y];
glint           = [params.glint.X, params.glint.Y];
viewDist        = params.viewDist;

%% Find the corner target furthest from the center
[~,I]           = max(nansum(abs(targets),2));
cT              = targets(I, :);
cP              = pupil(I, :);
cG              = glint(I, :);
%% Calculate the adjustment factor
rpc             = (sqrt((cT(1))^2 + (cT(2))^2 + viewDist^2) ...
    / sqrt((cT(1))^2 + (cT(2))^2)) * sqrt((cP(1) - cG(1))^2 + (cP(2) - cG(2))^2);
function [calMat] = calcCalMat(params)

%   Calculates the calibration matrix needed to determine gaze location 
%   from pupil and glint values in eye tracking data
%
%       The output transformation matrix is used to convert from pixels in 
%       an eye tracking video to mm on the display screen
%
%   Usage:
%       calMat = calcCalMat(params)
%
%   Inputs:
%       params.pupil.X      - vector of pupil X coordinates (pixels)
%       params.pupil.Y      - vector of pupil Y coordinates (pixels)
%       params.glint.X      - vector of glint X coordinates (pixels)
%       params.glint.Y      - vector of glint Y coordinates (pixels)
%       params.targets.X    - vector of target X coordinates (mm on screen)
%       params.targets.Y    - vector of target Y coordinates (mm on screen)
%       params.viewDist     - viewing distance (mm)
%       params.rpc          - adjustment parameter (see 'calcRpc')
%
%   Output:
%       calMat              - 4 x 4 transformation matrix
%
%   Written by Andrew S Bock Oct 2016

%% initialize the matrix
X = [...
    1 0 0 0 ...
    0 1 0 0 ...
    0 0 0 0 ...
    0 0 0 1 ...
    ];

%% exclude nan targets
params.pupil.X = params.pupil.X(~isnan(params.targets.X));
params.pupil.Y = params.pupil.Y(~isnan(params.targets.X));
params.glint.X = params.glint.X(~isnan(params.targets.X));
params.glint.Y = params.glint.Y(~isnan(params.targets.X));
params.targets.X = params.targets.X(~isnan(params.targets.X));
params.targets.Y = params.targets.Y (~isnan(params.targets.X));
%% Loop through calls to fminsearch, changing tolerance
for i=1:20
    options = optimset('Display','off','MaxFunEvals', 10000,...
        'MaxIter', 10000, 'TolX',10^(-i/2),'TolFun',10^(-i/2),'PlotFcns',[] );
    [X, f] = fminsearch(@(param) ...
        errfun(param,params.pupil,params.glint,params.targets,params.viewDist,params.rpc),...
        X, options);
    disp(['RSS error: ',num2str(f)])
end
% %% Loop through calls to fminunc, changing tolerance
% for i=1:30
%     options = optimoptions('fminunc', 'Algorithm','trust-region','Display','off', 'MaxFunctionEvaluations', 50000,...
%         'MaxIterations', 50000, 'OptimalityTolerance',10^(-i/2),'StepTolerance', 10^(-i/2),'PlotFcns',[] );
%     [X, f] = fminunc(@(param) ...
%         errfun(param,params.pupil,params.glint,params.targets,params.viewDist,params.rpc),...
%         X, options);
%     disp(['RSS error: ',num2str(f)])
% end
%% make the calibration matrix
calMat = [X(1:4); X(5:8); X(9:12); X(13:16)];

%% error function
function errtot = errfun(param, pupil, glint, targets, viewDist, rpc)

err = nan(1,length(targets(:,1)));
CalMatrix = [param(1:4); param(5:8); param(9:12); param(13:16)];

% minimize error for each target
for i = 1:length(targets.X)
    pX = pupil.X(i);
    pY = pupil.Y(i);
    gX = glint.X(i);
    gY = glint.Y(i);
    x = targets.X(i);
    y = targets.Y(i);
    z = viewDist;
    
    aXYZW = CalMatrix * [(pX-gX)/rpc; (pY-gY)/rpc; (1 - sqrt(((pX-gX)/rpc)^2 + ((pY-gY)/rpc)^2)); 1];
    oXYZW = [x; y; z; 1];
    
    errXYZ = (aXYZW(1:3)/aXYZW(4)) - (oXYZW(1:3)/oXYZW(4));
    err(i) = nansum(errXYZ.^2);
end
errtot = nansum(err.^2);
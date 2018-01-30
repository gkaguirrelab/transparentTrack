%% TEST_pupilProjectionModelInvertability
% Test if we can recover pupil azimuth, elevation, and radius from model
%
% Description:
%   transparentTrack implements a forward model of the appearance of the
%   pupil boundary in the image plane. This model begins with the azimuth
%   and elevation of the eye (in head-centered, fixed coordinates) and the
%   radius of the pupil. Using the parameters of a sceneGeometry structure
%   variable, the function `pupilProjection_fwd` returns the parameters of
%   an ellipse (in transparent form) in the image. A separate routine
%   (`pupilProjection_iv`) takes the parameters of an ellipse on the image
%   plane and (given a sceneGeometry structure) returns the eye parameters
%   of azimuth, elevation, and pupil radius.
%
%   Here we test, for a range of eye azimuths and elevations, how closely
%   we can recover the these values after passing through the forward and
%   inverse models.


%% Obtain a sceneGeometry structure
% If called with empty variables for the pupilData input and sceneGeometry
% output files, the estimateSceneGeometry routine returns a sceneGeometry
% file that contains the default parameters.
sceneGeometry = estimateSceneGeometry([],[]);

%% Obtain the ray tracing functions
rayTraceFuncs = assembleRayTraceFuncs( sceneGeometry );

%% Define some variables
pupilRadiusMM = 2;
eyeParams = [];
reconstructedEyeParams = [];
eyeParamErrors = [];

%% Loop over aimuths and elevations
% The range of values used here corresponds to the biological limits of the
% rotation of the eye horizontally and vertically.
thisTorsion = 0;
for thisAzimuth = -35:5:35
    for thisElevation = -25:5:25
        
        % Assemble the eyeParams variable
        eyeParams=[eyeParams; thisAzimuth,thisElevation,thisTorsion,pupilRadiusMM];
        
        % Forward projection from eyeParams to image ellipse
        tic
        pupilEllipseOnImagePlane = pupilProjection_fwd(eyeParams(end,:), sceneGeometry, rayTraceFuncs);
        toc
        
        % Inverse projection from image ellipse to eyeParams
        tic
        tmp = pupilProjection_inv(pupilEllipseOnImagePlane, sceneGeometry, rayTraceFuncs,'eyeParamsLB',[-40,-35,0,0.5],'eyeParamsUB',[40,35,0,4]);
        reconstructedEyeParams = [reconstructedEyeParams; tmp];
        toc
        
        % Calculate and save the error
        eyeParamErrors = [eyeParamErrors; eyeParams(end,:)-reconstructedEyeParams(end,:)];
    end
end

%% Report the errors
fprintf('The largest azimuth error is %f degrees.\n',max(eyeParamErrors(:,1)));
fprintf('The largest elevation error is %f degrees.\n',max(eyeParamErrors(:,2)));
fprintf('The largest radius error is %f millimeters.\n',max(eyeParamErrors(:,4)));
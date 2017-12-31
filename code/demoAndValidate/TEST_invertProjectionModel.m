%% UNIT TEST OF FORWARD AND INVERSE PUPIL PROJECTION MODEL

% Construct a sceneGeometry structure
sceneGeometry.intrinsicCameraMatrix = [772.5483 0 320; 0 772.5483 240; 0 0 1];
sceneGeometry.extrinsicRotationMatrix = [1 0 0; 0 -1 0; 0 0 -1];
sceneGeometry.eyeRadius = 11.29;
sceneGeometry.extrinsicTranslationVector = [0; 0; 50];
sceneGeometry.constraintTolerance = 0.02;

azimuths=[];
elevations=[];
azimuthErrors=[];
elevationErrors=[];
radiusErrors=[];

for thisAzimuth = -35:5:35
    for thisElevation = -25:5:25
        azimuths=[azimuths thisAzimuth];
        elevations=[elevations thisElevation];
            eyeParams=[thisAzimuth,thisElevation,2];
            pupilEllipseOnImagePlane = pupilProjection_fwd(eyeParams, sceneGeometry);
            reconstructedEyeParams = pupilProjection_inv(pupilEllipseOnImagePlane, sceneGeometry);
            eyeParamError = eyeParams-reconstructedEyeParams;
        azimuthErrors=[azimuthErrors eyeParamError(1)];
        elevationErrors=[elevationErrors eyeParamError(2)];
        radiusErrors=[radiusErrors eyeParamError(3)];
    end
end

fprintf('The largest azimuth error is: %f \n',max(azimuthErrors));
fprintf('The largest elevation error is: %f \n',max(elevationErrors));
fprintf('The largest radius error is: %f \n',max(radiusErrors));

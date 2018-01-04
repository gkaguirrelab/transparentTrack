%% UNIT TEST OF FORWARD AND INVERSE PUPIL PROJECTION MODEL

% Construct a sceneGeometry structure
sceneGeometry = estimateSceneGeometry([],[]);

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

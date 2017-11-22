%% UNIT TEST OF FORWARD AND INVERSE PUPIL PROJECTION MODELS

tolerance = 1e-9;

eyeCenter = [320 240 1500];
eyeRadius = 150;

rotationMax = 89;

projectionModels = {'pseudoPerspective' 'orthogonal' };

for models = 1:length(projectionModels)
    fprintf(['testing ' projectionModels{models} ' model \n']);
    for pupilAzimuth = -rotationMax:1:rotationMax
        for pupilElevation = -rotationMax:rotationMax:rotationMax
            reconstructedTransparentEllipse = pupilProjection_fwd(pupilAzimuth, pupilElevation, nan, eyeCenter, eyeRadius, projectionModels{models});
            [reconstructedPupilAzi, reconstructedPupilEle, ~] = pupilProjection_inv(reconstructedTransparentEllipse, eyeCenter, eyeRadius, projectionModels{models});
            if abs(reconstructedPupilAzi-pupilAzimuth) > tolerance || abs(reconstructedPupilEle-pupilElevation) > tolerance
                fprintf('Failed inversion check for azimuth %d, elevation %d \n',pupilAzimuth,pupilElevation);
                fprintf('   reconstruced azimuth %0.3f, reconstructed elevation %0.3f \n',reconstructedPupilAzi,reconstructedPupilEle);                
            end
        end
    end
end
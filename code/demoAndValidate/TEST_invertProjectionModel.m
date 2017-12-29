%% UNIT TEST OF FORWARD AND INVERSE PUPIL PROJECTION MODELS

tolerance = 1e-9;

eyeCenter = [320 240 1500];
eyeRadius = 150;
pupilArea = 1000;
projectionModels = {'pseudoPerspective' 'orthogonal' };

% %% Conduct some test
% 
% % Can we invert the forward projection? We test if the forward projection
% % ellipsed can be modeled with the sceneGeometry within a 0.1% error on
% % matching ellipse theta, eccentrivity, and error.
% for ii=1:length(pupilCenterAzimuths)
%     [eyeParams, bestMatchEllipseOnImagePlane, constraintViolations(ii)] = ...
%         pupilProjection_inv(projectedEllipsesOnImagePlane(ii,:), sceneGeometry, 'constraintTolerance', 0.001);
%     eyeParamError(ii,:)=eyeParams-[pupilCenterAzimuths(ii), pupilCenterElevations(ii), pupilRadii(ii)];
%     ellipseParamError(ii,:)=bestMatchEllipseOnImagePlane-projectedEllipsesOnImagePlane(ii,:);
% end
% fprintf('There were %d / %d constraint violations in the inversion of the forward model \n',sum(constraintViolations),length(pupilCenterAzimuths));


% Test if we can recover the azimuth and elevation of the eye after
% projecting the pupil to the image plane
allPass = true;
thetas=[];
rotationMax = 89;
for models = 1:length(projectionModels)
    fprintf(['testing ' projectionModels{models} ' model for eye --> ellipse --> eye \n']);
    for pupilAzimuth = -rotationMax:1:rotationMax
        for pupilElevation = -rotationMax:rotationMax:rotationMax
            reconstructedTransparentEllipse = pupilProjection_fwd(pupilAzimuth, pupilElevation, pupilArea, eyeCenter, eyeRadius, projectionModels{models});
            thetas=[thetas reconstructedTransparentEllipse(5)];
            [reconstructedPupilAzi, reconstructedPupilEle, reconstructedPupilArea] = pupilProjection_inv(reconstructedTransparentEllipse, eyeCenter, eyeRadius, projectionModels{models});
            if abs(reconstructedPupilAzi-pupilAzimuth) > tolerance || abs(reconstructedPupilEle-pupilElevation) > tolerance || abs(reconstructedPupilArea-pupilArea) > tolerance
                fprintf('Failed inversion check for azimuth %d, elevation %d \n',pupilAzimuth,pupilElevation);
                fprintf('   reconstruced azimuth %0.3f, reconstructed elevation %0.3f, reconstructed pupil area %0.3f \n',reconstructedPupilAzi,reconstructedPupilEle,reconstructedPupilArea);
                allPass=false;
            end
        end
    end
end

if allPass
    fprintf(['All inversion tests passed\n']);
end
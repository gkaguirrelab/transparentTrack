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
% createSceneGeometry returns a default sceneGeometry structure
sceneGeometry = createSceneGeometry();

%% Obtain the ray tracing functions
rayTraceFuncs = assembleRayTraceFuncs( sceneGeometry );

%% Define some variables
pupilRadiusMM = 2;
eyePoses = [];
reconstructedEyePoses = [];
eyePoseErrors = [];
centerErrors =[];
shapeErrors=[];
areaErrors=[];

%% Loop over aimuths and elevations
% The range of values used here corresponds to the biological limits of the
% rotation of the eye horizontally and vertically.
thisTorsion = 0;
for thisAzimuth = -35:5:35
    for thisElevation = -25:5:25
        
        % Assemble the eyePoses variable
        eyePoses=[eyePoses; thisAzimuth,thisElevation,thisTorsion,pupilRadiusMM];

        % Forward projection from eyePoses to image ellipse
        pupilEllipseOnImagePlane = pupilProjection_fwd(eyePoses(end,:), sceneGeometry, rayTraceFuncs);
        
        % Inverse projection from image ellipse to eyePoses. Note that we
        % must constrain at least one of the eye rotations, as the search
        % is otherwise unconstrained. We constrain torsion to be zero,
        % following Listing's Law.
        tic
        [inverseEyePose, bestMatchEllipseOnImagePlane, centerError, shapeError, areaError, exitFlag] = pupilProjection_inv(pupilEllipseOnImagePlane, sceneGeometry, rayTraceFuncs,'eyePoseLB',[-40,-35,0,0.5],'eyePoseUB',[40,35,0,4]);
        % If the exitFlag is 2, we may be in a local minimum. Repeat the 
        % search supplying the initial solution as x0.
        if exitFlag == 2
        [inverseEyePose, bestMatchEllipseOnImagePlane, centerError, shapeError, areaError, exitFlag] = pupilProjection_inv(pupilEllipseOnImagePlane, sceneGeometry, rayTraceFuncs,'eyePoseLB',[-40,-35,0,0.5],'eyePoseUB',[40,35,0,4],'x0',inverseEyePose);
        end
        toc
        reconstructedEyePoses = [reconstructedEyePoses; inverseEyePose];
        centerErrors=[centerErrors; centerError];
        shapeErrors=[shapeErrors; shapeError];
        areaErrors=[areaError; areaError];
        
        % Calculate and save the error
        eyePoseErrors = [eyePoseErrors; eyePoses(end,:)-reconstructedEyePoses(end,:)];
        
        if any(abs(eyePoseErrors(end,:))>0.05)
            foo = 1;
        end
    end
end

%% Report the errors
fprintf('The largest azimuth error is %f degrees.\n',max(abs(eyePoseErrors(:,1))));
fprintf('The largest elevation error is %f degrees.\n',max(abs(eyePoseErrors(:,2))));
fprintf('The largest radius error is %f millimeters.\n',max(abs(eyePoseErrors(:,4))));

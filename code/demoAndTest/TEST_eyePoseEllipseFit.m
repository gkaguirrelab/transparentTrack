%% TEST_eyePoseEllipseFit
% Test if we can recover pupil azimuth, elevation, and radius
%
% Description:
%   We attempt to recover the pose of the eye given a sceneGeometry and the
%   boundary points of an ellipse on the image plane


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
RMSEs =[];

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
        
        % Obtain boundary points for this ellipse
        [ Xp, Yp ] = ellipsePerimeterPoints( pupilEllipseOnImagePlane );
        
        % Inverse projection from image ellipse to eyePoses. Note that we
        % must constrain at least one of the eye rotations, as the search
        % is otherwise unconstrained. We constrain torsion to be zero,
        % following Listing's Law.
        tic
        [inverseEyePose, RMSE] = eyePoseEllipseFit(Xp, Yp, sceneGeometry, rayTraceFuncs,'eyePoseLB',[-40,-35,0,0.5],'eyePoseUB',[40,35,0,4]);
        toc
        
        reconstructedEyePoses = [reconstructedEyePoses; inverseEyePose];
        RMSEs = [RMSEs, RMSE];
        
        % Calculate and save the error
        eyePoseErrors = [eyePoseErrors; eyePoses(end,:)-reconstructedEyePoses(end,:)];
        
    end
end

%% Report the errors
fprintf('The largest azimuth error is %f degrees.\n',max(abs(eyePoseErrors(:,1))));
fprintf('The largest elevation error is %f degrees.\n',max(abs(eyePoseErrors(:,2))));
fprintf('The largest radius error is %f millimeters.\n',max(abs(eyePoseErrors(:,4))));

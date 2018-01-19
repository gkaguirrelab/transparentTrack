function displaceEyeCoordThroughCornea( sceneGeometry, eyeRotation )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% Build a rayTrace function
% This function traces a ray arising on the pupil through the cornea within
% a plane that includes the optic axis of the eye and a "height" dimension

    eye = sceneGeometry.eye;
    syms theta
    syms pupilPointHeight
    coords = [eye.pupilCenter(1) pupilPointHeight];
    opticalSystem = [nan nan eye.aqueousRefractiveIndex; ...
                     eye.corneaBackSurfaceCenter(1) -eye.corneaBackSurfaceRadius eye.corneaRefractiveIndex; ...
                     eye.corneaFrontSurfaceCenter(1) -eye.corneaFrontSurfaceRadius 1.0];
    outputRay = rayTraceCenteredSphericalSurfaces(coords, theta, opticalSystem);
    unitRayFromPupilFunc = matlabFunction(outputRay);
    clear outputRay theta pupilPointHeight
    
    % Create a function that takes a fixed pupil height and a symbolic
    % theta
    
    syms theta_p1p2
    syms theta_p1p3
    pupilPointHeight_p2 = 2;
    pupilPointHeight_p3 = 0;
    outputRayEyeWorld_p1p2 = unitRayFromPupilFunc(pupilPointHeight_p2, theta_p1p2);
    outputRayEyeWorld_p1p3 = unitRayFromPupilFunc(pupilPointHeight_p3, theta_p1p3);

    % Adjust the p1 (optical axis) position of the p1p3 ray to have the
    % same initial p1 positio nas the p1p2 ray.
    slope =outputRayEyeWorld_p1p3(2,2)/(outputRayEyeWorld_p1p3(2,1)-outputRayEyeWorld_p1p3(1,1));
    zOffset=outputRayEyeWorld_p1p2(1,1)-outputRayEyeWorld_p1p3(1,1);
    outputRayEyeWorld_p1p3(:,1)=outputRayEyeWorld_p1p3(:,1)+zOffset;
    outputRayEyeWorld_p1p3(:,2)=outputRayEyeWorld_p1p3(:,2)+(zOffset*slope);
    
    outputRayHeadWorld = (eyeRotation*(outputRayEyeWorld-sceneGeometry.eye.rotationCenter)')'+sceneGeometry.eye.rotationCenter;

    % Re-arrange the head world coordinate frame to transform to the scene
% world coordinate frame
    outputRaySceneWorld = outputRayHeadWorld(:,[2 3 1]);

    % We reverse the direction of the Y axis so that positive elevation of the
    % eye corresponds to a movement of the pupil upward in the image
    outputRaySceneWorld(:,2) = outputRaySceneWorld(:,2)*(-1);
    
    
    


end


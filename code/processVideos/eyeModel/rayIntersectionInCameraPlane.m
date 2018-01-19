function [zCameraPlaneX,zCameraPlaneY] = rayIntersectionInCameraPlane( sceneGeometry, eyeRotation, corneaRayTraceFunc )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

   
    
    syms theta_p1p2
    syms theta_p1p3
    syms pupilPointHeight_p2;
    syms pupilPointHeight_p3;
    outputRayEyeWorld_p1p2 = corneaRayTraceFunc(pupilPointHeight_p2, theta_p1p2);
    outputRayEyeWorld_p1p3 = corneaRayTraceFunc(pupilPointHeight_p3, theta_p1p3);

    % Adjust the p1 (optical axis) position of the p1p3 ray to have the
    % same initial p1 positio nas the p1p2 ray.
    slope =outputRayEyeWorld_p1p3(2,2)/(outputRayEyeWorld_p1p3(2,1)-outputRayEyeWorld_p1p3(1,1));
    zOffset=outputRayEyeWorld_p1p2(1,1)-outputRayEyeWorld_p1p3(1,1);
    outputRayEyeWorld_p1p3(:,1)=outputRayEyeWorld_p1p3(:,1)+zOffset;
    outputRayEyeWorld_p1p3(:,2)=outputRayEyeWorld_p1p3(:,2)+(zOffset*slope);
    
    outputRayEyeWorld(1,:) = [outputRayEyeWorld_p1p2(1,1) outputRayEyeWorld_p1p2(1,2) outputRayEyeWorld_p1p3(1,2)] - sceneGeometry.eye.rotationCenter;
    outputRayEyeWorld(2,:) = [outputRayEyeWorld_p1p2(2,1) outputRayEyeWorld_p1p2(2,2) outputRayEyeWorld_p1p3(2,2)] - sceneGeometry.eye.rotationCenter;
        
    outputRayHeadWorld = (eyeRotation*(outputRayEyeWorld)')';
    
    outputRayHeadWorld(1,:)=outputRayHeadWorld(1,:)+sceneGeometry.eye.rotationCenter;
    outputRayHeadWorld(2,:)=outputRayHeadWorld(2,:)+sceneGeometry.eye.rotationCenter;

    % Re-arrange the head world coordinate frame to transform to the scene
    % world coordinate frame
    outputRaySceneWorld = outputRayHeadWorld(:,[2 3 1]);

    % We reverse the direction of the Y axis so that positive elevation of the
    % eye corresponds to a movement of the pupil upward in the image
    outputRaySceneWorld(:,2) = outputRaySceneWorld(:,2)*(-1);
    
    % Obtain an expression for X and Y distances between the nodal point of the camera in the sceneWorld plane and the
    % point at which the ray will strike the plane that contains the camera
    slope_xZ =(outputRaySceneWorld(2,1)-outputRaySceneWorld(1,1))/(outputRaySceneWorld(2,3)-outputRaySceneWorld(1,3));
    slope_yZ =(outputRaySceneWorld(2,2)-outputRaySceneWorld(1,2))/(outputRaySceneWorld(2,3)-outputRaySceneWorld(1,3));

    zCameraPlaneX = outputRaySceneWorld(1,1)+((sceneGeometry.extrinsicTranslationVector(3)-outputRaySceneWorld(1,3))*slope_xZ);
    zCameraPlaneY = outputRaySceneWorld(1,2)+((sceneGeometry.extrinsicTranslationVector(3)-outputRaySceneWorld(1,3))*slope_yZ);
    
%    double(subs(zCameraPlaneX,theta_p1p2,pi/8))
%    double(subs(zCameraPlaneY,[ theta_p1p2, theta_p1p3],[pi/8,0.01]))
   

end


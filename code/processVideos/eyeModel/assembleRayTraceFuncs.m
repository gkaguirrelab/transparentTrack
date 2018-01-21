function [rayTraceFuncs] = assembleRayTraceFuncs( sceneGeometry )

% 2D ray tracing through the cornea
rayTraceFuncs.cornea = cornea(sceneGeometry);

% 2D distance of ray intersection on camera plane from camera node
[rayTraceFuncs.cameraNodeDistanceError2D.p1p2, rayTraceFuncs.cameraNodeDistanceError2D.p1p3] = ...
    cameraNodeDistanceError2D(rayTraceFuncs.cornea);

% 3D distance of ray intersection on camera plane from camera node
rayTraceFuncs.cameraNodeDistanceError3D = ...
    cameraNodeDistanceError3D(rayTraceFuncs.cornea);

% Virtual image ray
rayTraceFuncs.virtualImageRay = virtualImageRay(rayTraceFuncs.cornea);


end


%% LOCAL FUNCTIONS


function corneaFunc = cornea(sceneGeometry)
% 2D ray tracing through the cornea
% Assumes the camera resides in air (index of refraction is 1.0)
eye = sceneGeometry.eye;
syms z h theta
opticalSystem = [nan nan eye.aqueousRefractiveIndex; ...
    eye.corneaBackSurfaceCenter(1) -eye.corneaBackSurfaceRadius eye.corneaRefractiveIndex; ...
    eye.corneaFrontSurfaceCenter(1) -eye.corneaFrontSurfaceRadius 1.0];
outputRay2D = rayTraceCenteredSphericalSurfaces([z h], theta, opticalSystem);
corneaFunc = matlabFunction(outputRay2D);
end


function [p1p2Func, p1p3Func] = cameraNodeDistanceError2D(corneaFunc)
%% 2D Distance of the ray intersection from the camera node
syms p1 p2 p3
syms theta_p1p2 theta_p1p3

% The corneaFunc takes input as h, theta, z
outputRayEyeWorld2D_p1p2 = corneaFunc(p2,theta_p1p2,p1);
outputRayEyeWorld2D_p1p3 = corneaFunc(p3,theta_p1p3,p1);

% Add the third, constant dimension for the output rays
outputRayEyeWorld_p1p2=[outputRayEyeWorld2D_p1p2(1,1) outputRayEyeWorld2D_p1p2(1,2) p3;...
    outputRayEyeWorld2D_p1p2(2,1) outputRayEyeWorld2D_p1p2(2,2) p3];
outputRayEyeWorld_p1p3=[outputRayEyeWorld2D_p1p3(1,1) p2 outputRayEyeWorld2D_p1p3(1,2);...
    outputRayEyeWorld2D_p1p3(2,1) p2 outputRayEyeWorld2D_p1p3(2,2)];

% prepare to rotate the outputRay into the sceneWorld coordinates
syms cameraTranslationX cameraTranslationY cameraTranslationZ rotationCenterDepth
syms eyeAzimuthRads eyeElevationRads
eyeTorsionRads=0;

R3 = [cos(eyeAzimuthRads) -sin(eyeAzimuthRads) 0; sin(eyeAzimuthRads) cos(eyeAzimuthRads) 0; 0 0 1];
R2 = [cos(eyeElevationRads) 0 sin(eyeElevationRads); 0 1 0; -sin(eyeElevationRads) 0 cos(eyeElevationRads)];
R1 = [1 0 0; 0 cos(eyeTorsionRads) -sin(eyeTorsionRads); 0 sin(eyeTorsionRads) cos(eyeTorsionRads)];

eyeRotation = R1*R2*R3;

% Shift the eyeWorld ray to the rotational center of the eye,
% rotate for this eye pose, undo the centering
outputRayHeadWorld_p1p2(1,:)=outputRayEyeWorld_p1p2(1,:)-rotationCenterDepth;
outputRayHeadWorld_p1p2(2,:)=outputRayEyeWorld_p1p2(2,:)-rotationCenterDepth;
outputRayHeadWorld_p1p2 = (eyeRotation*(outputRayHeadWorld_p1p2)')';
outputRayHeadWorld_p1p2(1,:)=outputRayHeadWorld_p1p2(1,:)+rotationCenterDepth;
outputRayHeadWorld_p1p2(2,:)=outputRayHeadWorld_p1p2(2,:)+rotationCenterDepth;

outputRayHeadWorld_p1p3(1,:)=outputRayEyeWorld_p1p3(1,:)-rotationCenterDepth;
outputRayHeadWorld_p1p3(2,:)=outputRayEyeWorld_p1p3(2,:)-rotationCenterDepth;
outputRayHeadWorld_p1p3 = (eyeRotation*(outputRayHeadWorld_p1p3)')';
outputRayHeadWorld_p1p3(1,:)=outputRayHeadWorld_p1p3(1,:)+rotationCenterDepth;
outputRayHeadWorld_p1p3(2,:)=outputRayHeadWorld_p1p3(2,:)+rotationCenterDepth;

% Re-arrange the head world coordinate frame to transform to the scene
% world coordinate frame
outputRaySceneWorld_p1p2 = outputRayHeadWorld_p1p2(:,[2 3 1]);
outputRaySceneWorld_p1p3 = outputRayHeadWorld_p1p3(:,[2 3 1]);

% We reverse the direction of the Y axis so that positive elevation of the
% eye corresponds to a movement of the pupil upward in the image
outputRaySceneWorld_p1p2(:,2) = outputRaySceneWorld_p1p2(:,2)*(-1);
outputRaySceneWorld_p1p3(:,2) = outputRaySceneWorld_p1p3(:,2)*(-1);

% Obtain an expression for X and Y distances between the nodal point of the camera in the sceneWorld plane and the
% point at which the ray will strike the plane that contains the camera
slope_xZ =(outputRaySceneWorld_p1p2(2,1)-outputRaySceneWorld_p1p2(1,1))/(outputRaySceneWorld_p1p2(2,3)-outputRaySceneWorld_p1p2(1,3));
slope_yZ =(outputRaySceneWorld_p1p2(2,2)-outputRaySceneWorld_p1p2(1,2))/(outputRaySceneWorld_p1p2(2,3)-outputRaySceneWorld_p1p2(1,3));
cameraPlaneX = outputRaySceneWorld_p1p2(1,1)+((cameraTranslationZ-outputRaySceneWorld_p1p2(1,3))*slope_xZ);
cameraPlaneY = outputRaySceneWorld_p1p2(1,2)+((cameraTranslationZ-outputRaySceneWorld_p1p2(1,3))*slope_yZ);
p1p2Func = matlabFunction(sqrt(...
    (cameraTranslationX-cameraPlaneX)^2 + ...
    (cameraTranslationY-cameraPlaneY)^2 ...
    ));

slope_xZ =(outputRaySceneWorld_p1p3(2,1)-outputRaySceneWorld_p1p3(1,1))/(outputRaySceneWorld_p1p3(2,3)-outputRaySceneWorld_p1p3(1,3));
slope_yZ =(outputRaySceneWorld_p1p3(2,2)-outputRaySceneWorld_p1p3(1,2))/(outputRaySceneWorld_p1p3(2,3)-outputRaySceneWorld_p1p3(1,3));
cameraPlaneX = outputRaySceneWorld_p1p3(1,1)+((cameraTranslationZ-outputRaySceneWorld_p1p3(1,3))*slope_xZ);
cameraPlaneY = outputRaySceneWorld_p1p3(1,2)+((cameraTranslationZ-outputRaySceneWorld_p1p3(1,3))*slope_yZ);
p1p3Func = matlabFunction(sqrt(...
    (cameraTranslationX-cameraPlaneX)^2 + ...
    (cameraTranslationY-cameraPlaneY)^2 ...
    ));

end


function p1p2p3Func = cameraNodeDistanceError3D(corneaFunc)
%% 3D Distance of the ray intersection from the camera node
syms p1 p2 p3
syms theta_p1p2 theta_p1p3

% The corneaFunc takes input as h, theta, z
outputRayEyeWorld2D_p1p2 = corneaFunc(p2,theta_p1p2,p1);
outputRayEyeWorld2D_p1p3 = corneaFunc(p3,theta_p1p3,p1);


% Create a 3D output ray system. Shift the p1p3 ray to have the same
% initial p1 value as the p1p2 vector
slope =(outputRayEyeWorld2D_p1p2(2,2)-outputRayEyeWorld2D_p1p2(1,2))/(outputRayEyeWorld2D_p1p2(2,1)-outputRayEyeWorld2D_p1p2(1,1));
zOffset=outputRayEyeWorld2D_p1p2(1,1)-p1;
outputRayEyeWorld2D_p1p2(:,1)=outputRayEyeWorld2D_p1p2(:,1)-zOffset;
outputRayEyeWorld2D_p1p2(:,2)=outputRayEyeWorld2D_p1p2(:,2)-(zOffset*slope);

slope =(outputRayEyeWorld2D_p1p3(2,2)-outputRayEyeWorld2D_p1p3(1,2))/(outputRayEyeWorld2D_p1p3(2,1)-outputRayEyeWorld2D_p1p3(1,1));
zOffset=outputRayEyeWorld2D_p1p3(1,1)-p1;
outputRayEyeWorld2D_p1p3(:,1)=outputRayEyeWorld2D_p1p3(:,1)-zOffset;
outputRayEyeWorld2D_p1p3(:,2)=outputRayEyeWorld2D_p1p3(:,2)-(zOffset*slope);

% Combine into a single, 3D ray
outputRayEyeWorld3D=[outputRayEyeWorld2D_p1p2(1,1) outputRayEyeWorld2D_p1p2(1,2) outputRayEyeWorld2D_p1p3(1,2);...
    outputRayEyeWorld2D_p1p2(2,1) outputRayEyeWorld2D_p1p2(2,2) outputRayEyeWorld2D_p1p3(2,2)];

% prepare to rotate the outputRay into the sceneWorld coordinates
syms cameraTranslationX cameraTranslationY cameraTranslationZ rotationCenterDepth
syms eyeAzimuthRads eyeElevationRads
eyeTorsionRads=0;

R3 = [cos(eyeAzimuthRads) -sin(eyeAzimuthRads) 0; sin(eyeAzimuthRads) cos(eyeAzimuthRads) 0; 0 0 1];
R2 = [cos(eyeElevationRads) 0 sin(eyeElevationRads); 0 1 0; -sin(eyeElevationRads) 0 cos(eyeElevationRads)];
R1 = [1 0 0; 0 cos(eyeTorsionRads) -sin(eyeTorsionRads); 0 sin(eyeTorsionRads) cos(eyeTorsionRads)];

eyeRotation = R1*R2*R3;

% Shift the eyeWorld ray to the rotational center of the eye,
% rotate for this eye pose, undo the centering
outputRayHeadWorld3D(1,:)=outputRayEyeWorld3D(1,:)-rotationCenterDepth;
outputRayHeadWorld3D(2,:)=outputRayEyeWorld3D(2,:)-rotationCenterDepth;
outputRayHeadWorld3D = (eyeRotation*(outputRayHeadWorld3D)')';
outputRayHeadWorld3D(1,:)=outputRayHeadWorld3D(1,:)+rotationCenterDepth;
outputRayHeadWorld3D(2,:)=outputRayHeadWorld3D(2,:)+rotationCenterDepth;


% Re-arrange the head world coordinate frame to transform to the scene
% world coordinate frame
outputRaySceneWorld3D = outputRayHeadWorld3D(:,[2 3 1]);

% We reverse the direction of the Y axis so that positive elevation of the
% eye corresponds to a movement of the pupil upward in the image
outputRaySceneWorld3D(:,2) = outputRaySceneWorld3D(:,2)*(-1);

% Obtain an expression for X and Y distances between the nodal point of the camera in the sceneWorld plane and the
% point at which the ray will strike the plane that contains the camera
slope_xZ =(outputRaySceneWorld3D(2,1)-outputRaySceneWorld3D(1,1))/(outputRaySceneWorld3D(2,3)-outputRaySceneWorld3D(1,3));
slope_yZ =(outputRaySceneWorld3D(2,2)-outputRaySceneWorld3D(1,2))/(outputRaySceneWorld3D(2,3)-outputRaySceneWorld3D(1,3));
cameraPlaneX = outputRaySceneWorld3D(1,1)+((cameraTranslationZ-outputRaySceneWorld3D(1,3))*slope_xZ);
cameraPlaneY = outputRaySceneWorld3D(1,2)+((cameraTranslationZ-outputRaySceneWorld3D(1,3))*slope_yZ);
p1p2p3Func = matlabFunction(sqrt(...
    (cameraTranslationX-cameraPlaneX)^2 + ...
    (cameraTranslationY-cameraPlaneY)^2 ...
    ));

end



function [virtualImageRayFunc] = virtualImageRay(corneaFunc)

syms p1 p2 p3
syms theta_p1p2 theta_p1p3

% Assemble the virtual image ray
outputRayEyeWorld_p1p2 = corneaFunc(p2, theta_p1p2, p1);
outputRayEyeWorld_p1p3 = corneaFunc(p3, theta_p1p3, p1);

% Adjust the p1 (optical axis) position of the rays to have the
% their initial position at the same p1
slope =(outputRayEyeWorld_p1p2(2,2)-outputRayEyeWorld_p1p2(1,2))/(outputRayEyeWorld_p1p2(2,1)-outputRayEyeWorld_p1p2(1,1));
zOffset=outputRayEyeWorld_p1p2(1,1)-p1;
outputRayEyeWorld_p1p2(:,1)=outputRayEyeWorld_p1p2(:,1)-zOffset;
outputRayEyeWorld_p1p2(:,2)=outputRayEyeWorld_p1p2(:,2)-(zOffset*slope);

slope =(outputRayEyeWorld_p1p3(2,2)-outputRayEyeWorld_p1p3(1,2))/(outputRayEyeWorld_p1p3(2,1)-outputRayEyeWorld_p1p3(1,1));
zOffset=outputRayEyeWorld_p1p3(1,1)-p1;
outputRayEyeWorld_p1p3(:,1)=outputRayEyeWorld_p1p3(:,1)-zOffset;
outputRayEyeWorld_p1p3(:,2)=outputRayEyeWorld_p1p3(:,2)-(zOffset*slope);

% Combine the two dimensions into a single, 3D ray
outputRayEyeWorld(1,:) = [outputRayEyeWorld_p1p2(1,1) outputRayEyeWorld_p1p2(1,2) outputRayEyeWorld_p1p3(1,2)];
outputRayEyeWorld(2,:) = [outputRayEyeWorld_p1p2(2,1) outputRayEyeWorld_p1p2(2,2) outputRayEyeWorld_p1p3(2,2)];

virtualImageRayFunc = matlabFunction(outputRayEyeWorld);
end

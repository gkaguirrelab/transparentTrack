function corneaRayTraceFunc = createCorneaRayTraceFunction( sceneGeometry )

    % Create the 2D ray tracing function through the cornea. Assumes the
    % output medium in which the camera resides is air (index of
    % refraction is 1.0)

    eye = sceneGeometry.eye;
    syms z h theta
    opticalSystem = [nan nan eye.aqueousRefractiveIndex; ...
                     eye.corneaBackSurfaceCenter(1) -eye.corneaBackSurfaceRadius eye.corneaRefractiveIndex; ...
                     eye.corneaFrontSurfaceCenter(1) -eye.corneaFrontSurfaceRadius 1.0];
    outputRay2D = rayTraceCenteredSphericalSurfaces([z h], theta, opticalSystem);
    corneaRayTraceFunc2D = matlabFunction(outputRay2D);
    
    % Inspect this variable to confirm that the order of input variables
    % for the corneaRayTraceFunc2D is height, theta, z
    varsCRTF2D = symvar(outputRay2D);
    
    clear z h theta outputRay2D opticalSystem eye

    % Assemble the 3D ray tracing function

    syms p1 p2 p3
    syms theta_p1p2 theta_p1p3
    outputRayEyeWorld_p1p2 = corneaRayTraceFunc2D(p2, theta_p1p2, p1);
    outputRayEyeWorld_p1p3 = corneaRayTraceFunc2D(p3, theta_p1p3, p1);

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

    % Combine the two dimensions into a single, 3D function
    outputRayEyeWorld(1,:) = [outputRayEyeWorld_p1p2(1,1) outputRayEyeWorld_p1p2(1,2) outputRayEyeWorld_p1p3(1,2)];
    outputRayEyeWorld(2,:) = [outputRayEyeWorld_p1p2(2,1) outputRayEyeWorld_p1p2(2,2) outputRayEyeWorld_p1p3(2,2)];

    corneaRayTraceFunc = matlabFunction(outputRayEyeWorld);
    
end
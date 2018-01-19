function [ cornealRayTraceFunc ] = createCornealRayTraceFunction( sceneGeometry )

% Assumes the output medium in which the camera resides is air (index of
% refraction is 1.0)

    eye = sceneGeometry.eye;
    syms theta
    syms pupilPointHeight
    coords = [eye.pupilCenter(1) pupilPointHeight];
    opticalSystem = [nan nan eye.aqueousRefractiveIndex; ...
                     eye.corneaBackSurfaceCenter(1) -eye.corneaBackSurfaceRadius eye.corneaRefractiveIndex; ...
                     eye.corneaFrontSurfaceCenter(1) -eye.corneaFrontSurfaceRadius 1.0];
    outputRay = rayTraceCenteredSphericalSurfaces(coords, theta, opticalSystem);
    cornealRayTraceFunc = matlabFunction(outputRay);

end


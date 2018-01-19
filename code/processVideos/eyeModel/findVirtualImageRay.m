function findVirtualImageRay( eyeWorldPoints, sceneGeometry, eyeRotation )


[ corneaRayTraceFunc ] = createCorneaRayTraceFunction( sceneGeometry );
[zCameraPlaneX,zCameraPlaneY] = rayIntersectionInCameraPlane( sceneGeometry, eyeRotation, corneaRayTraceFunc );

syms pupilPointHeight_p2 theta_p1p2 theta_p1p3
 k=subs(zCameraPlaneX,[ pupilPointHeight_p2, theta_p1p2],[2, theta_p1p2]);
eqn = subs(k,theta_p1p2,theta_p1p2) ==0;
solve(eqn,theta_p1p2)

end


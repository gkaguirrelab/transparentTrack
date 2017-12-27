function [distanceErrorVec, eyeParams] = crazyFunction( ellipses, sceneGeometry )

% Identify the center of projection. We use this information to constrain
% the eye azimuths and elevations that we will allow to fit a given ellipse
projectionMatrix = ...
    sceneGeometry.intrinsicCameraMatrix * ...
    [sceneGeometry.extrinsicRotationMatrix, ...
    sceneGeometry.extrinsicTranslationVector];

CoP = projectionMatrix*[0 0 0 1]';
CoP(1:2)=CoP(1:2)./CoP(3);
CoP=CoP(1:2);


options = optimoptions(@fmincon,...
    'Display','off',...
    'Algorithm','interior-point',...
    'Diagnostics','off');


ellipseAreaError = @(e1,e2) abs(e1(3)-e2(3));
ellipseDistanceError = @(e1,e2) sqrt((e1(1)-e2(1))^2 + (e1(2)-e2(2))^2);

for ii = 1:size(ellipses,1)

    myNonLinFunc = @(p) constrainEccentrictyAndTheta( pupilProjection_fwd(p(1), p(2), p(3), sceneGeometry), ellipses(ii,:) );
    myErrorFunc = @(p) ellipseAreaError( pupilProjection_fwd(p(1), p(2), p(3), sceneGeometry), ellipses(ii,:) );
    if ellipses(ii,1) < CoP(1)
        AziLB = -35; AziUB = 0;
    else
        AziLB = 0; AziUB = 35;
    end
    if ellipses(ii,2) > CoP(2)
        EleLB = -25; EleUB = 0;
    else
        EleLB = 0; EleUB = 25;
    end
        
    p = fmincon(myErrorFunc,[0 0 2],[],[],[],[],[AziLB EleLB 1],[AziUB EleUB 4],myNonLinFunc, options);
    distanceErrorVec(ii) = ellipseDistanceError(pupilProjection_fwd(p(1), p(2), p(3), sceneGeometry), ellipses(ii,:));
    eyeParams(ii,:)=p;
end

end
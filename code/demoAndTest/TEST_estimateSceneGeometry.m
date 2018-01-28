%% TEST_estimateSceneGeometry
% Examine the ability of the routines to estimate an unknown scene geometry
%
% Description:

%

close all

% Obtain the default sceneGeometry
% - no lens distortion
% - camera distance to 100 mm
% - center of rotation at the corneal apex
%
% By having the eye rotate around the corneal apex, we match the conditions
% of these prior papers in which the camera was maintained at a constant
% distance from the corneal apex.
veridicalSceneGeometry = estimateSceneGeometry([],[],'eyeLaterality','Right');

% Assemble the ray tracing functions
rayTraceFuncs = assembleRayTraceFuncs( veridicalSceneGeometry );

% Create a set of ellipses from the veridial geometry
idx=1;
for azi=-15:15:15
    for ele=-15:15:15
        pupilRadius = 2*randn()/2;
        eyeParams=[azi, ele, 0, pupilRadius];
        pupilData.initial.ellipses.values(idx,:) = pupilProjection_fwd(eyeParams, veridicalSceneGeometry, rayTraceFuncs);
        pupilData.initial.ellipses.RMSE(idx) = 1;
        idx=idx+1;
    end
end

estimatedSceneGeometry = estimateSceneGeometry(pupilData,'','useParallel',true,'ellipseArrayList',1:1:idx-1);



%% TEST_pupilProjection_inv
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
%   inverse models. We also examine the effect of not performing ray
%   tracing in the inverse projection.


%% Obtain a sceneGeometry structure
% createSceneGeometry returns a default sceneGeometry structure
sceneGeometry = createSceneGeometry();

%% Obtain the ray tracing functions
rayTraceFuncs = assembleRayTraceFuncs( sceneGeometry );

%% Define some variables
pupilRadiusMM = 2;
eyePoses = [];
pupilEllipseAreas = zeros(15,11);
eyePoseErrorsWithRayTrace = zeros(15,11,4);
eyePoseErrorsWithoutRayTrace = zeros(15,11,4);


%% Loop over aimuths and elevations
% The range of values used here corresponds to the biological limits of the
% rotation of the eye horizontally and vertically.
for aziIdx = 1:15
    for eleIdx = 1:11
        
        thisAzimuth = (aziIdx-8)*5;
        thisElevation = (eleIdx-6)*5;
        thisTorsion = 0;
        
        % Assemble the eyePoses variable
        eyePose=[eyePoses; thisAzimuth,thisElevation,thisTorsion,pupilRadiusMM];
        
        % Forward projection from eyePoses to image ellipse
        pupilEllipseOnImagePlane = pupilProjection_fwd(eyePose, sceneGeometry, rayTraceFuncs);
        pupilEllipseAreas(aziIdx,eleIdx) = pupilEllipseOnImagePlane(3);
        
        % Inverse projection from image ellipse to eyePoses. Note that we
        % must constrain at least one of the eye rotations, as the search
        % is otherwise underdetermined. We constrain torsion to be zero,
        % following Listing's Law.
        tic
        [inverseEyePose, ~, ~, ~, ~, exitFlag] = ...
            pupilProjection_inv(pupilEllipseOnImagePlane, sceneGeometry, rayTraceFuncs,'eyePoseLB',[-40,-35,0,0.5],'eyePoseUB',[40,35,0,4]);
        % If the exitFlag is 2, we may be in a local minimum. Repeat the
        % search, supplying the initial solution as x0.
        if exitFlag == 2
            x0 = inverseEyePose + [1e-3 1e-3 0 1e-3];
            [inverseEyePose] = ...
                pupilProjection_inv(pupilEllipseOnImagePlane, sceneGeometry, rayTraceFuncs,'eyePoseLB',[-40,-35,0,0.5],'eyePoseUB',[40,35,0,4],'x0',x0);
        end
        toc
        % Save the error
        eyePoseErrorsWithRayTrace(aziIdx,eleIdx,:) = eyePose-inverseEyePose;
        
        
        % Perform the search again, this time without ray tracing.
        tic
        [inverseEyePose, ~, ~, ~, ~, exitFlag] = ...
            pupilProjection_inv(pupilEllipseOnImagePlane, sceneGeometry, [],'eyePoseLB',[-40,-35,0,0.5],'eyePoseUB',[40,35,0,4]);
        % If the exitFlag is 2, we may be in a local minimum. Repeat the
        % search, supplying the initial solution as x0.
        if exitFlag == 2
            x0 = inverseEyePose + [1e-3 1e-3 0 1e-3];
            [inverseEyePose] = ...
                pupilProjection_inv(pupilEllipseOnImagePlane, sceneGeometry, [],'eyePoseLB',[-40,-35,0,0.5],'eyePoseUB',[40,35,0,4],'x0',x0);
        end
        toc
        % Save the error
        eyePoseErrorsWithoutRayTrace(aziIdx,eleIdx,:) = eyePose-inverseEyePose;
    end
end

%% Report the errors
fprintf('The largest azimuth error is %f degrees.\n',max(max(abs(eyePoseErrorsWithRayTrace(:,:,1)))));
fprintf('The largest elevation error is %f degrees.\n',max(max(abs(eyePoseErrorsWithRayTrace(:,:,2)))));
fprintf('The largest proportion radius error is %f.\n',1-max(max(abs(eyePoseErrorsWithRayTrace(:,:,4)-pupilRadiusMM)./pupilRadiusMM)));

%% Create some figures
idxToPlot = [1,2,4];
plotRange = [-5 5; -5 5; -0.05 0.05];
titleStrings = {'azimuth error','elevation error','-proportion pupil radius error'};

figure
for panel = 1:3
    subplot(3,1,panel)
    if panel == 3
        image = squeeze(eyePoseErrorsWithoutRayTrace(:,:,idxToPlot(panel)))+pupilRadiusMM;
        % Assert that the pupil size at
        image = -((image - image(8,6))./image(8,6))';
    else
        image = squeeze(eyePoseErrorsWithoutRayTrace(:,:,idxToPlot(panel)))';
    end
    [nr,nc] = size(image);
    pcolor([image nan(nr,1); nan(1,nc+1)]);
    caxis(plotRange(panel,:));
    shading flat;
    axis equal
    % Set the axis backgroud to dark gray
    set(gcf,'Color',[1 1 1]); set(gca,'Color',[.75 .75 .75]); set(gcf,'InvertHardCopy','off');
    colorbar;
    title(titleStrings{panel});
    xlabel('veridical azimuth [deg]')
    ylabel('veridical elevation [deg]')
    xticks((1:1:size(image,2))+.5);
    xticklabels(-35:5:35);
    xtickangle(90);
    yticks((1:1:size(image,1))+.5);
    yticklabels(-25:5:25);
    xlim([1 size(image,2)+1]);
    ylim([1 size(image,1)+1]);
end

figure
image = -((pupilEllipseAreas-pupilEllipseAreas(8,6))./pupilEllipseAreas(8,6))';
[nr,nc] = size(image);
pcolor([image nan(nr,1); nan(1,nc+1)]);
caxis([-0.125 0.125]);
shading flat;
axis equal
% Set the axis backgroud to dark gray
set(gcf,'Color',[1 1 1]); set(gca,'Color',[.75 .75 .75]); set(gcf,'InvertHardCopy','off');
colorbar;
title('-proportion pupil area error');
xlabel('veridical azimuth [deg]')
ylabel('veridical elevation [deg]')
xticks((1:1:size(image,2))+.5);
xticklabels(-35:5:35);
xtickangle(90);
yticks((1:1:size(image,1))+.5);
yticklabels(-25:5:25);
xlim([1 size(image,2)+1]);
ylim([1 size(image,1)+1]);

% TEST_sceneGeometryConstraint

% This script generates a video composed of ellipses on the image plane
% that are pseudoPerspective projections of a model pupil. The video is
% then run through the pupil analysis pipeline, including derivation of the
% scene geometry and sceneConstrained ellipse fitting.

%% Clean up
clearvars
close all

% Reset the random number generator so that we have reproducible results
rng default

% define scene dimension in pixels
videoSizeX = 640;
videoSizeY = 480;

% Setup the video save location and path params
sandboxDir = '~/Desktop/sceneGeometryConstraintTEST';
pathParams.dataOutputDirFull = fullfile(sandboxDir);
pathParams.dataSourceDirFull = fullfile(sandboxDir);
pathParams.runName = 'synthetic';
videoName = fullfile(sandboxDir, [pathParams.runName '_gray.avi']);

% check or make a directory for output
if exist(sandboxDir,'dir')==0
    mkdir(sandboxDir);
end


%% construct scene geometry
% The Center of Rotation is aligned with the center of the image plane
% along the optical axis of the camera. The pupil rotates on a 125 pixel
% radius. The pupil has a fixed radius as sweeps the scene back and forth
% at different elevations.
eyeRadius =  332.0173;
pupilRadius = 59.1135;
sceneDistance = 805.7827;

% Set up the sceneGeometry
projectionModel = 'pseudoPerspective';
eyeCenter = [videoSizeX/2, videoSizeY/2, eyeRadius+sceneDistance];


% Define the range of azimuth and elevation steps
eleSteps = -20:5:20;
aziSweeps = -30:3:30;

allPupilAzi=[];
allPupilEle=[];
for ii = 1: length(eleSteps)
    allPupilEle = [allPupilEle eleSteps(ii)*ones(1,length(aziSweeps))];
    if mod(ii,2)
        allPupilAzi = [allPupilAzi aziSweeps];
    else
        allPupilAzi = [allPupilAzi fliplr(aziSweeps)];
    end
end

%% create the video
% the video will be created with a modified version of the
% pupilProjection_fwd function, as the Y axis for matlab plots is oriented
% upwards. Normally pupilProjection_fwd operates on video frames, for which
% the convention is that the Y axis points downwards.

nFrames = numel(allPupilAzi);
emptyFrame = ones(videoSizeY, videoSizeX)*0.85;
h = figure('visible','off');

% for each frame, create the projected ellipse and test the consistenpFitExplicit(2) of
% the returned theta values
for ii = 1:nFrames
    % derive ellipse points from Azimuth and Elevation values
    pupilAzi = allPupilAzi(ii); % in degrees
    pupilEle = allPupilEle(ii); % in degrees
    
    forwardProjectEllipseParams = pupilProjection_fwd_Yup(pupilAzi, pupilEle, pi*pupilRadius.^2, eyeCenter, eyeRadius, projectionModel);
    pFitExplicit = ellipse_transparent2ex(forwardProjectEllipseParams);
    
    % make the plot and save it as a frame
    tempImage=emptyFrame;
    fh=@(x,y) pFitExplicit(4)^2.*((x-pFitExplicit(1)).*cos(pFitExplicit(5))-(y-pFitExplicit(2)).*sin(pFitExplicit(5))).^2 + pFitExplicit(3)^2.*((x-pFitExplicit(1)).*sin(pFitExplicit(5))+(y-pFitExplicit(2)).*cos(pFitExplicit(5))).^2 - pFitExplicit(3)^2.*pFitExplicit(4)^2 < 0;
    [Y, X]=ind2sub(size(tempImage),1:1:numel(tempImage));
    tempImage(fh(X,Y))=0;
    
    % overlay a glint
    glintPositionX = forwardProjectEllipseParams(1) - (forwardProjectEllipseParams(1)-videoSizeX/2).*0.15;
    glintPositionY = forwardProjectEllipseParams(2) - (forwardProjectEllipseParams(2)-videoSizeY/2).*0.15;    
%      glintPositionX = videoSizeX/2 + 2 ;
%     glintPositionY = videoSizeY/2 + 2;  
    tempImage = insertShape(tempImage,'filledCircle',[glintPositionX, glintPositionY, 5],'Color','w','Opacity',1);
    
    % display the frame
    imshow(tempImage,'Border','tight');
    axis equal
    axis off
    ylim([0 videoSizeY]);
    xlim([0 videoSizeX]);
    truesize;
    
    % store the frame
    thisFrame(ii) = getframe(gca);
    hold off
    
end % loop over frames

% Create a color map for the indexed video
cmap = [linspace(0,1,256)' linspace(0,1,256)' linspace(0,1,256)'];
cmap(1,:)=[1 0 0];
cmap(2,:)=[0 1 0];
cmap(3,:)=[0 0 1];
cmap(4,:)=[1 1 0];
cmap(5,:)=[0 1 1];
cmap(6,:)=[1 0 1];

% write the video
writerObj = VideoWriter(videoName, 'Uncompressed AVI');
writerObj.FrameRate = 60;
open(writerObj);
for ii=1:nFrames
    indexedFrame = rgb2ind(thisFrame(ii).cdata,cmap, 'nodither');
    writeVideo(writerObj,indexedFrame);
end
close (writerObj);


%% Perform the analysis with one call
% runVideoPipeline( pathParams, ...
%     'verbosity', 'full', 'useParallel',false, 'catchErrors', false,...
%     'maskBox', [0.9 0.9], ...
%     'overwriteControlFile',true, ...
%     'ellipseTransparentUB',[videoSizeX,videoSizeY, 20000, 1.0, pi],...
%     'eyeRadius',eyeRadius, 'cameraDistanceInPixels',sceneDistance, ...
%     'sceneGeometryLB',[0, 0, sceneDistance+eyeRadius, 25],'sceneGeometryUB',[640, 480, sceneDistance+eyeRadius, 500],...
%     'skipStageByNumber',1);

%% run short version of analysis
perimeterFileName = fullfile(sandboxDir, 'synthetic_perimeter.mat');
pupilFileName = fullfile(sandboxDir, 'syntheticPerimeter_pupil.mat');
sceneGeometryFileName = fullfile(sandboxDir, 'syntheticPerimeter_sceneGeometry.mat');
sceneDiagnosticPlotFileName = fullfile(sandboxDir, 'syntheticPerimeter_sceneDiagnosticPlot.pdf');
finalFitVideoName = fullfile(sandboxDir, 'syntheticPerimeter_finalFit.avi');


findPupilPerimeter(videoName,perimeterFileName,'verbosity','full','maskBox', [0.9 0.9]);
fitPupilPerimeter(perimeterFileName, pupilFileName,'verbosity','full','ellipseTransparentLB',[],'ellipseTransparentUB',[],'nSplits',0);
sceneGeometry = estimateSceneGeometry(pupilFileName, sceneGeometryFileName,'sceneDiagnosticPlotFileName', sceneDiagnosticPlotFileName,'sceneDiagnosticPlotSizeXY', [videoSizeX videoSizeY], ...
    'projectionModel','pseudoPerspective', ...
    'eyeRadius',eyeRadius, 'cameraDistanceInPixels',sceneDistance, ...
    'sceneGeometryLB',[0, 0, sceneDistance, 25],'sceneGeometryUB',[640, 480, sceneDistance, 500], ...
    'verbosity','full');
fitPupilPerimeter(perimeterFileName,pupilFileName,'sceneGeometryFileName',sceneGeometryFileName,'ellipseTransparentLB',[0, 0, 300, 0, 0],'ellipseTransparentUB',[videoSizeX,videoSizeY,20000,1.0, pi],'verbosity','full');
pupilData = smoothPupilArea(perimeterFileName, pupilFileName, sceneGeometryFileName,'verbosity','full');
makeFitVideo(videoName, finalFitVideoName, 'pupilFileName',pupilFileName,'sceneGeometryFileName',sceneGeometryFileName,'perimeterFileName',perimeterFileName,'perimeterColor','r','whichFieldToPlot','ellipseParamsAreaSmoothed_mean','verbosity','full')

%% Verify that the scene geometry allows for the correct reconstruction of the eye position
% note on the pupil area reconstruction: matlab uses different pixel
% indexing conventions when creating plots and analizing videos. For
% plotting, the origin of the XY cartesian plane is conventionally put at
% (0,0). For image processing (hence for the whole pupil tracking pipeline)
% the origin of the pixel is put at (0.5,0.5). This causes the pupil radius
% to be calculated half a pixel shorter during the analysis. We account for
% that definining a pupilRadiusOnImage, which is half a pixel shorter than
% pupilRadius we used to construct the TEST video.

% load(fullfile(pathParams.dataOutputDirFull,[pathParams.runName '_pupil.mat']));
% load(fullfile(pathParams.dataOutputDirFull,[pathParams.runName '_sceneGeometry.mat']));

ellipses = pupilData.ellipseParamsAreaSmoothed_mean;
pupilRadiusOnImage = pupilRadius - 0.5;

for ii = 1:nFrames
    [reconstructedPupilAzi(ii), reconstructedPupilEle(ii), reconstructedPupilArea(ii)] = ...
        pupilProjection_inv(ellipses(ii,:),  [sceneGeometry.eyeCenter.X sceneGeometry.eyeCenter.Y, sceneGeometry.eyeCenter.Z], sceneGeometry.eyeRadius, sceneGeometry.meta.projectionModel);
end

% Report how closely we have reconstructed the actual scene geometry
fprintf('Veridical scene geometry - eye center: [%0.1f, %0.1f, %0.1f], eye radius: %0.1f \n',eyeCenter(1), eyeCenter(2), eyeCenter(3), eyeRadius);
fprintf('Estimated scene geometry - eye center: [%0.1f, %0.1f, %0.1f], eye radius: %0.1f \n',sceneGeometry.eyeCenter.X, sceneGeometry.eyeCenter.Y, sceneGeometry.eyeCenter.Z, sceneGeometry.eyeRadius);

% plot real Azi,  Ele, and pupil area vs reconstructed ones
figure
subplot(1,3,1)
plot(allPupilAzi,reconstructedPupilAzi, '.r')
rl = refline(1,0);
rl.Color = 'k';
xlabel('Ground Truth Pupil Azimuth in degrees')
ylabel('Reconstructed Pupil Azimuth in degrees')
ylim([min([min(allPupilAzi) min(reconstructedPupilAzi)]) max([max(allPupilAzi) max(reconstructedPupilAzi)])]);
xlim([min([min(allPupilAzi) min(reconstructedPupilAzi)]) max([max(allPupilAzi) max(reconstructedPupilAzi)])]);
axis square

subplot(1,3,2)
plot(-allPupilEle,reconstructedPupilEle, '.r')
rl = refline(1,0);
rl.Color = 'k';
xlabel('Ground Truth Pupil Elevation in degrees')
ylabel('Reconstructed Pupil Elevation in degrees')
ylim([min([min(allPupilAzi) min(reconstructedPupilAzi)]) max([max(allPupilAzi) max(reconstructedPupilAzi)])]);
xlim([min([min(allPupilAzi) min(reconstructedPupilAzi)]) max([max(allPupilAzi) max(reconstructedPupilAzi)])]);
axis square

subplot(1,3,3)
plot(1:1:nFrames,ellipses(:,3), '.b')
hold on
plot(1:1:nFrames,reconstructedPupilArea, 'xr')
line(1:1:nFrames, (pi*(pupilRadiusOnImage)^2)*ones(1,nFrames),'Color', 'black');
xlabel('frame')
ylabel('reconstructed pupil area')
axis square





%% LOCAL FUNCTION

function reconstructedTransparentEllipse = pupilProjection_fwd_Yup(pupilAzi, pupilEle, pupilArea, eyeCenter, eyeRadius, projectionModel)
% modified version of pupilProjection_fwd, to take into account that the Y
% axis points UP when constructing the model (in the videos the Y axis
% points conventionally down)

%% main

% initiate reconstructedTransparentEllipse
reconstructedTransparentEllipse = nan(1,5);

% if we have a non-defined case, just return all nans
if any(isnan([pupilAzi pupilEle]))
    return
end

% calculate the pupilCenter3D
pupilCenter3D(1) = ( eyeRadius) *sind(pupilAzi)*cosd(pupilEle);
pupilCenter3D(2) = -(eyeRadius)*sind(pupilEle);
pupilCenter3D(3) = (eyeCenter(3) - eyeRadius)-(eyeRadius)*cosd(pupilAzi)*cosd(pupilEle);

% define ellipse center
switch projectionModel
    case 'orthogonal'
        % under orthogonal hypothesis the ellipse center in 2D is
        % coincident with the ellipse center in the plane of projection.
        reconstructedTransparentEllipse(1) = pupilCenter3D(1) + eyeCenter(1);
        reconstructedTransparentEllipse(2) = pupilCenter3D(2) + eyeCenter(2);
    case 'pseudoPerspective'
        % for the pseudoPerspective correction, we uniformly scale the
        % orthogonal projection according to the scene distance and the eye
        % radius.
        
        % get the perspective projection correction factor
        sceneDistance = eyeCenter(3) - eyeRadius;
        perspectiveCorrectionFactor = sceneDistance/(sceneDistance + pupilCenter3D(3));
        
        % apply perspective correction factor to the pupil center
        reconstructedTransparentEllipse(1) = (pupilCenter3D(1) * perspectiveCorrectionFactor) + eyeCenter(1);
        reconstructedTransparentEllipse(2) = (pupilCenter3D(2) * perspectiveCorrectionFactor) + eyeCenter(2);
end

% derive ellipse eccentricity
e = sqrt((sind(pupilEle))^2 - ((sind(pupilAzi))^2 * ((sind(pupilEle))^2 - 1)));
reconstructedTransparentEllipse(4) = e;

% derive ellipse theta
theta = nan;
if pupilAzi > 0  &&  pupilEle > 0
    theta = - asin(sind(pupilAzi)/e);
elseif pupilAzi > 0 && pupilEle < 0
    theta =  asin(sind(pupilAzi)/e);
elseif pupilAzi < 0  &&  pupilEle > 0
    theta = - asin(sind(pupilAzi)/e);
elseif pupilAzi < 0 && pupilEle < 0
    theta =  asin(sind(pupilAzi)/e);
elseif abs(pupilAzi) < 1e-12 && abs(pupilEle) < 1e-12
    theta = 0;
elseif abs(pupilAzi) < 1e-12 && abs(pupilEle) >= 1e-12
    theta = 0;
elseif abs(pupilEle) < 1e-12 && abs(pupilAzi) >= 1e-12
    theta = pi/2;
else
    % Couldn't constrain the theta. This shouldn't happen
    warning('For some reason the theta was unconstrained. Setting to nan');
end

% Keep the theta values between 0 and pi
if theta < 0
    theta = theta+pi;
end
if theta > pi
    theta = theta - pi;
end

reconstructedTransparentEllipse(5) = theta;

% area (if pupilArea available)
if ~isnan(pupilArea)
    switch projectionModel
        case 'orthogonal'
            % under orthogonal hypothesis, the semimajor axis of the
            % ellipse equals the pupil circle radius. If that is known, it
            % can be assigned and used later to determine the ellipse area.
            pupilRadius = sqrt(pupilArea/pi);
            semiMajorAxis = pupilRadius;
        case 'pseudoPerspective'
            % apply scaling factor to pupil area
            pupilRadius = sqrt(pupilArea/pi);
            semiMajorAxis = pupilRadius.*perspectiveCorrectionFactor;
    end
    semiMinorAxis = semiMajorAxis*sqrt(1-e^2);
    reconstructedTransparentEllipse(3) = pi*semiMajorAxis*semiMinorAxis;
end

end % function


function XYZ = applyRoll(xyz, rotationDeg)
    x = xyz(1); y = xyz(2); z=xyz(3);
	X = x;
    Y = y .* cosd(rotationDeg) - z .* sind(rotationDeg);
    Z = y .* sind(rotationDeg) + z .* cosd(rotationDeg);
    XYZ=[X Y Z];
end

function XYZ = applyPitch(xyz, rotationDeg)
    x = xyz(1); y = xyz(2); z=xyz(3);
	X = x .* cosd(rotationDeg) + z .* sind(rotationDeg);
    Y = y;
    Z = - x .* sind(rotationDeg) + z .* cosd(rotationDeg);
    XYZ=[X Y Z];
end

function XYZ = applyYaw(xyz, rotationDeg)
    x = xyz(1); y = xyz(2); z=xyz(3);
    X = x .* cosd(rotationDeg) - y .* sind(rotationDeg);
    Y = x .* sind(rotationDeg) + y .* cosd(rotationDeg);
    Z = z;
    XYZ=[X Y Z];
end


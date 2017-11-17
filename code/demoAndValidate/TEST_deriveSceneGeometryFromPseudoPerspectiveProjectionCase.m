% TEST_deriveSceneGeometryFromPseudoPerspectiveProjectionCase

% this script generates a video recording of a "pupil" moving in 3D space.
% The pupil is centered in the scene, has a fixed area and rotates around a
% fixed length radius. The rotations are defined by an Azimuth and
% an Elevation angle.
% The video is saved and later analyzed using the standard processing
% function up to reconstructing the scene geometry.
% Finally the pupil is projected back in 3D space using the ellipse fits
% and the reconstructed scene geometry, to verify that the reconstruction
% matches the ground truth movement pattern.

%% Clean up
clearvars
close all


%% Define sandbox dir
sandboxDir = '~/Desktop/sceneGeometryFromPseudoPerspectiveTEST';

% check or make a directory for output
    if exist(sandboxDir,'dir')==0
        mkdir(sandboxDir);
    end

%% UnitVideo settings
syntheticPerimVideoName = fullfile(sandboxDir, 'synthetic.avi');

% 1. define synthetic data length
nFrames = 200;

% 2. define scene dimension in pixels
videoSizeX = 640;
videoSizeY = 480;

% create the video writer object
writerObj = VideoWriter(syntheticPerimVideoName, 'Uncompressed AVI');
writerObj.FrameRate = 60;


% define eye movements
% define rotations in deg
allPupilAzi = [randn(1,nFrames)*15]; % in degrees
allPupilEle = [randn(1,nFrames)*15]; % in degrees


%% construct scene geometry
% we will construct a scene where the CoR sits right at the center, and the
% pupil rotates on a 150 pixels radius.

% define the eye sphere radius
eyeballR =  150+5;
eyeballCenter = [0 0 0];
eyeball = [eyeballCenter eyeballR];

% "Depth" of the pupil plane, with respect to the sphere radius. The
% deeper, the bigger the pupil.
planeDepth = 5;
intersectingSphere = [eyeballCenter eyeballR-planeDepth];

% rotation arm length
rotationArmLength = eyeballR - planeDepth;

% create scene plane
sceneDistance = 1200; % orthogonal distance from rotation arm
xMax= videoSizeX/2; % max x size of scene (for plotting purposes)
yMax = videoSizeY/2; % max y size of scene
scenePlane = createPlane([0 0 rotationArmLength+sceneDistance],[0 0 rotationArmLength+sceneDistance]);

% center of projection
centerOfProjection3D = projPointOnPlane(eyeballCenter, scenePlane);

% derive optical axis
% [centerOfProjectionX,centerOfProjectionY,centerOfProjectionZ] = sph2cart(azi0,ele0,rotationArmLength);
opticalAxis = createLine3d(eyeballCenter,centerOfProjection3D);


%% create the scene and save a frame in the video
% we save a video here, where the background is white and the pupil ellipse
% is black and filled.

emptyFrame = ones(videoSizeY, videoSizeX);

h = figure('visible','off');
% create scene, save a frame
for ii = 1:length(allPupilAzi)
    pupilAzi = allPupilAzi(ii); % in degrees
    pupilEle = allPupilEle(ii); % in degrees
    
    % first rotate along Y (azi)
    rotationX = createRotationOy(eyeballCenter, deg2rad(pupilAzi));
    pupilAxis = transformLine3d(opticalAxis, rotationX);
    
    % then along X (elevation)
    rotationY = createRotationOx(eyeballCenter, deg2rad(-pupilEle)); % - to preserve convention
    pupilAxis = transformLine3d(pupilAxis, rotationY);
    
    % intersect sphere and pupil axis
    pupilCenter3D = intersectLineSphere(pupilAxis, intersectingSphere);
    
    % create the pupil plane
    pupilPlane = createPlane(pupilCenter3D(2,:),pupilCenter3D(2,:)-eyeballCenter);
    
    % intersect pupil plane and sphere to obtain the pupil circle
    pupilInEye = intersectPlaneSphere(pupilPlane, eyeball);
    pupilRadius = pupilInEye(4); % saving it explicitely for reference
    
    %  get the 3d cartesian coordinates of the pupil perimeter points
    pupilPoints3d = getCirclePoints3d(pupilInEye);
    
    % do orthogonal projeciton on the scene plane
    pupilPoints2d = projPointOnPlane(pupilPoints3d,scenePlane);
    
    % add a little noise to the points
    pupilPoints2d = awgn(pupilPoints2d,.5);
    
    % calculate the perspective correction factor
    relativeDepth = rotationArmLength*(1-(cosd(pupilEle)*cosd(pupilAzi)));
    perspectiveCorrectionFactor = (sceneDistance/(sceneDistance + relativeDepth))
    
    % scale the pupil points according to the perspective correction factor
    pupilPoints2d = pupilPoints2d .* perspectiveCorrectionFactor;
                
    % make the plot and save it as a frame
    imshow(emptyFrame);
    hold on
    fill(pupilPoints2d(:,1)+xMax,pupilPoints2d(:,2)+yMax,'k')
    axis equal
    axis off
    ylim([0 videoSizeY]);
    xlim([0 videoSizeX]);
    truesize;
    hold off
    
    % get the frame
    thisFrame(ii) = getframe(gca);
    
end

% Create a color map for the indexed video
cmap = [linspace(0,1,256)' linspace(0,1,256)' linspace(0,1,256)'];
cmap(1,:)=[1 0 0];
cmap(2,:)=[0 1 0];
cmap(3,:)=[0 0 1];
cmap(4,:)=[1 1 0];
cmap(5,:)=[0 1 1];
cmap(6,:)=[1 0 1];

% write the video
open(writerObj);

for ii=1:nFrames
    indexedFrame = rgb2ind(thisFrame(ii).cdata,cmap, 'nodither');
    
    
    % get video size
    videoX = size(indexedFrame,2);
    videoY = size(indexedFrame,1);
    
    % write video
    writeVideo(writerObj,indexedFrame);
end

close (writerObj);

%% use this video in the video pipeline 
syntheticPerimFileName = fullfile(sandboxDir, 'synthetic_perimeter.mat');
pupilFileName = fullfile(sandboxDir, 'syntheticPerimeter_pupil.mat');
sceneGeometryFileName = fullfile(sandboxDir, 'syntheticPerimeter_sceneGeometry.mat');
sceneDiagnosticPlotFileName = fullfile(sandboxDir, 'syntheticPerimeter_sceneDiagnosticPlot.pdf');
finalFitVideoName = fullfile(sandboxDir, 'syntheticPerimeter_finalFit.avi');

findPupilPerimeter(syntheticPerimVideoName,syntheticPerimFileName,'verbosity','full');
pupilData = fitPupilPerimeter(syntheticPerimFileName, pupilFileName,'verbosity','full','ellipseTransparentLB',[0, 0, 300, 0, -0.5*pi],'ellipseTransparentUB',[videoX,videoY,20000,0.75, 0.5*pi],'nSplits',0);
sceneGeometry = estimateSceneGeometry(pupilFileName, sceneGeometryFileName,'sceneDiagnosticPlotFileName', sceneDiagnosticPlotFileName,'sceneDiagnosticPlotSizeXY', [videoX videoY], ...
    'projectionModel','pseudoPerspective','eyeRadius',rotationArmLength, 'cameraDistanceInPixels',sceneDistance);

%% Verify that the scene geometry allows for the correct reconstruction of the eye movements

ellipses = pupilData.ellipseParamsUnconstrained_mean;
eyeCenter = [sceneGeometry.eyeCenter.X sceneGeometry.eyeCenter.Y, sceneGeometry.eyeCenter.Z];
projectionModel = 'pseudoPerspective';

for ii = 1:nFrames
    [reconstructedPupilAzi(ii), reconstructedPupilEle(ii), reconstructedPupilArea(ii)] = pupilProjection_inv(ellipses(ii,:),  eyeCenter, rotationArmLength, projectionModel);
end

% plot real Azi and Ele vs reconstructed ones
figure
subplot(1,2,1)
plot(allPupilAzi,reconstructedPupilAzi, '.')
xlabel('Ground Truth Pupil Azimuth in degrees')
ylabel('Reconstructed Pupil Azimuth in degrees')

subplot(1,2,2)
plot(allPupilEle,reconstructedPupilEle, '.')
xlabel('Ground Truth Pupil Elevation in degrees')
ylabel('Reconstructed Pupil Elevation in degrees')
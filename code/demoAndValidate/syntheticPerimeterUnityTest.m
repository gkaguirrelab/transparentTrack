% script UNITY TEST 2


% make an artificial perimeter video, based on a "real eye"

sandboxDir = '~/Desktop/eyeTrackingDEMO';
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
allPupilAzi = [0 randn(1,nFrames-1)*15]; % in degrees
allPupilEle = [0 randn(1,nFrames-1)*15]; % in degrees


%% construct scene geometry
% we will construct a scene where the CoR sits right at the center, and the
% pupil rotates on a 250 pixels radius. The first ellipse is, by
% construction, projected as a true circle.

% define the eye sphere radius
eyeballR =  250+15;
eyeballCenter = [0 0 0];
eyeball = [eyeballCenter eyeballR];

% "Depth" of the pupil plane, with respect to the sphere radius. The
% deeper, the bigger the pupil.
planeDepth = 5;
intersectingSphere = [eyeballCenter eyeballR-planeDepth];


% rotation arm length
rotationArmLength = eyeballR - planeDepth;

% create scene plane
sceneDistance = 100; % orthogonal distance from rotation arm
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
%    pupilPoints2d = awgn(pupilPoints2d,3);
    
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

%% use this video in the pipeline
syntheticPerimFileName = fullfile(sandboxDir, 'synthetic_perimeter.mat');
pupilFileName = fullfile(sandboxDir, 'syntheticPerimeter_pupil.mat');
sceneGeometryFileName = fullfile(sandboxDir, 'syntheticPerimeter_sceneGeometry.mat');
sceneDiagnosticPlotFileName = fullfile(sandboxDir, 'syntheticPerimeter_sceneDiagnosticPlot.pdf');
finalFitVideoName = fullfile(sandboxDir, 'syntheticPerimeter_finalFit.avi');

findPupilPerimeter(syntheticPerimVideoName,syntheticPerimFileName,'verbosity','full');
fitPupilPerimeter(syntheticPerimFileName, pupilFileName,'verbosity','full','ellipseTransparentLB',[0, 0, 300, 0, -0.5*pi],'ellipseTransparentUB',[videoX,videoY,20000,0.75, 0.5*pi],'nSplits',0);
estimateSceneGeometry(pupilFileName, sceneGeometryFileName,'sceneDiagnosticPlotFileName', sceneDiagnosticPlotFileName,'sceneDiagnosticPlotSizeXY', [videoX videoY]);

% script UNITY TEST 2


% make an artificial perimeter video, based on a "real eye"

sandboxDir = '~/Desktop/eyeTrackingDEMO';
%% UnitVideo settings
syntheticPerimVideoName = fullfile(sandboxDir, 'syntheticPerim.avi');
syntheticPerimFileName = fullfile(sandboxDir, 'syntheticPerimeter.mat');
% 1. define synthetic data length
nFrames = 100;

% 2. define frame dimention in pixels
videoSizeX = 640;
videoSizeY = 480;

% create the video writer object
writerObj = VideoWriter(syntheticPerimVideoName, 'Uncompressed AVI');
writerObj.FrameRate = 60;



% define eye movements
% define rotations in deg
allPupilAzi = linspace(0,50,nFrames); % in degrees
allPupilEle = zeros(1,nFrames); % in degrees


%% construct scene geometry

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
xMax= 320; % max x size of scene (for plotting purposes)
yMax = 240; % max y size of scene
scenePlane = createPlane([0 0 rotationArmLength+sceneDistance],[0 0 rotationArmLength+sceneDistance]);

% center of projection
centerOfProjection3D = projPointOnPlane(eyeballCenter, scenePlane);


% derive optical axis
% [centerOfProjectionX,centerOfProjectionY,centerOfProjectionZ] = sph2cart(azi0,ele0,rotationArmLength);
opticalAxis = createLine3d(eyeballCenter,centerOfProjection3D);


%% create the scene and save a frame in the video

% initialize variable to hold the perimeter data
% perimeter_data = zeros(videoSizeY,videoSizeX,nFrames,'uint8');

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
    
    % add noise to the points
    pupilPoints2d = awgn(pupilPoints2d,2);
    
    % make the plot and save it as a frame    
    plot(pupilPoints2d(:,1),pupilPoints2d(:,2),'w')
    set(gca,'Color','k','PlotBoxAspectRatio',[4 3 1])
%     set(gca,'Units','pixels', 'Position',[0 0 videoSizeX videoSizeY])
    xlim([-xMax xMax])
    ylim([-yMax yMax])
    
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
    indexedFrame = imresize(indexedFrame,[480 640]);
    writeVideo(writerObj,indexedFrame);
    binP= im2bw(indexedFrame);
    perimFrame = im2uint8(binP);
    perimeter.data(:,:,ii) = perimFrame;
end

close (writerObj);
save(syntheticPerimFileName,'perimeter');

%% use this perimeter video in the pipeline

pupilFileName = fullfile(sandboxDir, 'syntheticPerimeter_pupil.mat');
sceneGeometryFileName = fullfile(sandboxDir, 'syntheticPerimeter_sceneGeometry.mat');
sceneDiagnosticPlotFileName = fullfile(sandboxDir, 'syntheticPerimeter_sceneDiagnosticPlot.pdf');
finalFitVideoName = fullfile(sandboxDir, 'syntheticPerimeter_finalFit.avi');

fitPupilPerimeter(syntheticPerimFileName, pupilFileName,'verbosity','full','ellipseTransparentLB',[0, 0, 500, 0, -0.5*pi],'ellipseTransparentUB',[640,480,20000,0.75, 0.5*pi],'nSplits',0);
estimateSceneGeometry(pupilFileName, sceneGeometryFileName,'sceneDiagnosticPlotFileName', sceneDiagnosticPlotFileName);

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
    
    % make the plot and save it as a frame    
    plot(pupilPoints2d(:,1),pupilPoints2d(:,2),'w')
    set(gca,'Color','k','PlotBoxAspectRatio',[4 3 1])
%     set(gca,'Units','pixels', 'Position',[0 0 videoSizeX videoSizeY])
    xlim([-xMax xMax])
    ylim([-yMax yMax])
    
    thisFrame = getframe(gca);
    perimeter_data(:,:,:,ii) = thisFrame.cdata;
end

% Create a color map for the indexed video
cmap = [linspace(0,1,256)' linspace(0,1,256)' linspace(0,1,256)'];
cmap(1,:)=[1 0 0];
cmap(2,:)=[0 1 0];
cmap(3,:)=[0 0 1];
cmap(4,:)=[1 1 0];
cmap(5,:)=[0 1 1];
cmap(6,:)=[1 0 1];

% save out the perimeter file
perimeter.data = perimeter_data;
save(syntheticPerimFileName,'perimeter');

% write the video
open(writerObj);

for ii=1:nFrames
    indexedFrame = rgb2ind(squeeze(perimeter_data(:,:,:,ii)), cmap, 'nodither');
    writeVideo(writerObj,indexedFrame);
end

close (writerObj);



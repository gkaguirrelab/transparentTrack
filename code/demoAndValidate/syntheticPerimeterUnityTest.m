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
allPupilAzi = linspace(0,60,nFrames); % in degrees
allPupilEle = zeros(1,nFrames); % in degrees


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
% by default, matlab will plot the scene with a 1:1 pixel ratio.

% To try and understand if pixel ratio influences the appearence of the
% pupil on the scene, we can comment in/out the PlotBoxAspectRatio line,
% that would set the plot pixel ratio to 4:3. This is the pixel ratio
% of the camera we use. We don't know for sure this is the pixel ratio of the
% video when it is digitalized by the V.TOP, but it is likely.


% select pixel ratio
pixelRatio = '4:3';
% pixelRatio = '1:1';

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
    pupilPoints2d = awgn(pupilPoints2d,3);
    
    % make the plot and save it as a frame
    fill(pupilPoints2d(:,1),pupilPoints2d(:,2),'k')
    
    % format the axis with a 4:3 ratio (default 1:1)
    switch pixelRatio
        case '4:3'
            set(gca,'PlotBoxAspectRatio',[4 3 1])
        case '1:1'
            % do nothing
    end
    xlim([-xMax xMax])
    ylim([-yMax yMax])
    
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
    
    switch pixelRatio
        case '1:1'
            % doubling the size of this image makes the size dimensions even
            % numbers (use for 1:1 case)
            indexedFrame = imresize(indexedFrame,2);
            
        case '4:3'
            % if the aspect ratio is corrected to 4:3, the video can be scaled to
            % match the LiveTrack+V.TOP video resolution.
            indexedFrame = imresize(indexedFrame,[videoSizeY videoSizeX]);
    end
    
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

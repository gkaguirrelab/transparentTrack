% DEMO find distances and COP

%% housekeeping
close all
clear all
%clc

%% load pupil data
pupilFileName = '~/Desktop/eyeTrackingDEMO/TOME_processing/session2_spatialStimuli/TOME_3020/050517/EyeTracking/GazeCal01_pupil.mat';

%% find sceneGeometry
sceneGeometryFileName = '~/Desktop/eyeTrackingDEMO/TOME_processing/session2_spatialStimuli/TOME_3020/050517/EyeTracking/GazeCal01_sceneGeometry.mat';
sceneGeometry = estimateSceneGeometry (pupilFileName,sceneGeometryFileName,'verbosity','full');

%% given the scene geometry and ellipses centers, find the eccentricity and theta to be used as fitting constraints

% load ellises from the pupilfile
load(pupilFileName)
ellipses = pupilData.pInitialFitTransparent;
% initialize eccentricity and theta
eccentricity = nan(length(ellipses),1);
theta = nan(length(ellipses),1);
for ii = 1:length(ellipses)
    if any(isnan(ellipses(ii)))
        continue
    else
        [eccentricity(ii,:), theta(ii)] = constrainEllipseBySceneGeometry (ellipses(ii,:),sceneGeometry);
    end
end

figure
plot(ellipses(:,4),eccentricity,'.')
xlabel('transparentEllipses eccentricities')
ylabel('constrained eccentricities')

figure
plot(ellipses(:,4))
hold on
plot(eccentricity)
title('Timeseries for eccentricity')
legend('transparentEllipses eccentricities','constrained eccentricities')

figure
plot(ellipses(:,5),theta,'.')
xlabel('transparentEllipses tilt angles')
ylabel('constrained tilt angles')

figure
plot(sin(ellipses(:,5)))
hold on
plot(sin(theta))
title('timeseries for the sin of the thetas') % for easier comparison
legend('transparentEllipses tilt angles','constrained tilt angles')
ylabel('sin(theta)')

%% use the derived eccentricity and thetas to reconstruct the ellipses
% we take the initalFit ellipses, maintain the X Y area values and
% substitue the constraint eccentricity and theta. Then, we plot this on
% top of the original video to see how the constrained fit looks like.
% (while it won't be perfect, we hope that it is not too far off)
%
constrainedEllipses = [ellipses(:,1:3) eccentricity theta];
grayVideoName = '~/Desktop/eyeTrackingDEMO/TOME_processing/session2_spatialStimuli/TOME_3020/050517/EyeTracking/GazeCal01_gray.avi';
constrainedEllipsesVideoName  = '~/Desktop/eyeTrackingDEMO/TOME_processing/session2_spatialStimuli/TOME_3020/050517/EyeTracking/GazeCal01_constrainedEllipses.avi';

% read video file into memory
videoInObj = VideoReader(grayVideoName);

% get number of frames
nFrames = floor(videoInObj.Duration*videoInObj.FrameRate);
nFrames = 1000;

% get video dimensions
videoSizeX = videoInObj.Width;
videoSizeY = videoInObj.Height;
% read the video into memory
for ii = 1:nFrames
    thisFrame = readFrame(videoInObj);
    sourceVideo(:,:,ii) = rgb2gray (thisFrame);
end
% close the video object
clear videoInObj

% prepare the outputVideo
outputVideo=zeros(videoSizeY,videoSizeX,3,nFrames,'uint8');

for ii = 1:nFrames
    
    % Create a figure
    frameFig = figure( 'Visible', 'off');
    
    % show the initial frame
    imshow(squeeze(sourceVideo(:,:,ii)), 'Border', 'tight');
    hold on
    
    if sum(isnan(constrainedEllipses(ii,:)))==0
        % build ellipse impicit equation
        pFitImplicit = ellipse_ex2im(ellipse_transparent2ex(constrainedEllipses(ii,:)));
        fh=@(x,y) pFitImplicit(1).*x.^2 +pFitImplicit(2).*x.*y +pFitImplicit(3).*y.^2 +pFitImplicit(4).*x +pFitImplicit(5).*y +pFitImplicit(6);
        % superimpose the ellipse using fimplicit or ezplot
        if exist('fimplicit','file')==2
            fimplicit(fh,[1, videoSizeX, 1, videoSizeY],'Color', 'red','LineWidth',1);
            set(gca,'position',[0 0 1 1],'units','normalized')
            axis off;
        else
            plotHandle=ezplot(fh,[1, videoSizeX, 1, videoSizeY]);
            set(plotHandle, 'Color', 'red')
            set(plotHandle,'LineWidth',1);
        end
    end
    
    % Save the frame and close the figure
    tmp=getframe(frameFig);
    outputVideo(:,:,:,ii)=tmp.cdata;
    close(frameFig);
    
end

% Create a color map
cmap = [linspace(0,1,256)' linspace(0,1,256)' linspace(0,1,256)'];
cmap(1,:)=[1 0 0];
cmap(2,:)=[0 1 0];
cmap(3,:)=[0 0 1];
cmap(4,:)=[1 1 0];
cmap(5,:)=[0 1 1];
cmap(6,:)=[1 0 1];

% write the outputVideo to file
videoOutObj = VideoWriter(constrainedEllipsesVideoName,'Indexed AVI');
videoOutObj.FrameRate = 60;
videoOutObj.Colormap = cmap;
open(videoOutObj);

% loop through the frames and save them
for ii=1:nFrames
    indexedFrame = rgb2ind(squeeze(outputVideo(:,:,:,ii)), cmap, 'nodither');
    writeVideo(videoOutObj,indexedFrame);
end
% close the videoObj
clear videoOutObj


%% derive pupil params in 3D
% pupil params in 3D are azimut and elevation (with respect to the center
% of projection on the scene) and pupil area
pupil3DFileName = '~/Desktop/eyeTrackingDEMO/TOME_processing/session2_spatialStimuli/TOME_3020/050517/EyeTracking/GazeCal01_pupil3D.mat';

for ii = 1:length(constrainedEllipses)
    centerOfProjection = [sceneGeometry.eyeballCenter.X sceneGeometry.eyeballCenter.Y]; % because of orthogonal projection
    [reconstructedPupilAzi(ii), reconstructedPupilEle(ii), reconstructedPupilRadius(ii)] = pupilProjection_inv(constrainedEllipses(ii,:),centerOfProjection);
end

pupil3D.azi = reconstructedPupilAzi;
pupil3D.ele = reconstructedPupilEle;
pupil3D.radius= reconstructedPupilRadius;
pupil3D.centerOfProjection = centerOfProjection;
save(pupil3DFileName,'pupil3D')
%% make a pupil movement animation reproducing the 3D movements of the eye

% first attempt: pupil of fixed size (just to see if this works)

synteticFixedSizePupilVideoName = '~/Desktop/eyeTrackingDEMO/TOME_processing/session2_spatialStimuli/TOME_3020/050517/EyeTracking/GazeCal01_syntFixedSizePupil.avi';



% 1.initial geometry

% define the eyeball sphere radius
eyeballR =  250;
eyeballCenter = [pupil3D.centerOfProjection(1) pupil3D.centerOfProjection(2) 0];
eyeball = [eyeballCenter eyeballR];

% "Depth" of the pupil plane, with respect to the sphere radius. The
% deeper, the bigger the pupil.
planeDepth = 15;
intersectingSphere = [eyeballCenter eyeballR-planeDepth];


% "Depth" of the intersecting plane, with respect to the sphere radius. The
% deeper, the bigger the pupil.
planeDepth = 5;

% rotation arm length
rotationArmLength = eyeballR - planeDepth;

% create scene plane
sceneDistance = 5000; % orthogonal distance from rotation arm
xMax= 640; % max x size of scene (for plotting purposes)
yMax = 480; % max y size of scene
scenePlane = createPlane([0 0 rotationArmLength+sceneDistance],[0 0 rotationArmLength+sceneDistance]);

% center of projection
centerOfProjection3D = projPointOnPlane(eyeballCenter, scenePlane);

% derive optical axis
% [centerOfProjectionX,centerOfProjectionY,centerOfProjectionZ] = sph2cart(azi0,ele0,rotationArmLength);
opticalAxis = createLine3d(eyeballCenter,centerOfProjection3D);


for ii = 1:length(pupil3D.azi)
    % define rotations in deg
    pupilAzi = pupil3D.azi(ii); % in degrees
    pupilEle = pupil3D.ele(ii); % in degrees
    
    if ~isnan(pupilAzi)
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
    end
    
    % Create a figure
    frameFig = figure( 'Visible', 'off');
    
    
    % add the eyeball
    drawSphere(eyeball,'FaceAlpha',0.5)
    hold on
    
    % add optical axis
    drawLine3d(opticalAxis,'k')
    hold on
    
    if ~isnan(pupilAzi)
        % add the pupil on the eye
        drawCircle3d(pupilInEye)
        hold on
        
        % add pupil axis
        drawLine3d(pupilAxis,'b')
        hold on
    end
    
    
    % add labels
    xlabel('X')
    ylabel('Y')
    zlabel('Z')
    
    view([-35 60])
    title( 'Reconstruction of 3D pupil movements (fixed size pupil)')
    axis equal
    
    
    
    % Save the frame and close the figure
    tmp=getframe(frameFig);
    outputVideo(:,:,:,ii)=tmp.cdata;
    close(frameFig);
end

% Create a color map
cmap = [linspace(0,1,256)' linspace(0,1,256)' linspace(0,1,256)'];
cmap(1,:)=[1 0 0];
cmap(2,:)=[0 1 0];
cmap(3,:)=[0 0 1];
cmap(4,:)=[1 1 0];
cmap(5,:)=[0 1 1];
cmap(6,:)=[1 0 1];

% write the outputVideo to file
videoOutObj = VideoWriter(synteticFixedSizePupilVideoName,'Indexed AVI');
videoOutObj.FrameRate = 60;
videoOutObj.Colormap = cmap;
open(videoOutObj);

% loop through the frames and save them
for ii=1:nFrames
    indexedFrame = rgb2ind(squeeze(outputVideo(:,:,:,ii)), cmap, 'nodither');
    writeVideo(videoOutObj,indexedFrame);
end
% close the videoObj
clear videoOutObj

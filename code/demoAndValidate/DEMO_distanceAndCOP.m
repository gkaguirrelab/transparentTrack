% DEMO find distances and COP

%% housekeeping
close all
clear all
%clc

%% load pupil data
pupilFileName = '~/Desktop/eyeTrackingDEMO/TOME_processing/session2_spatialStimuli/TOME_3020/050517/EyeTracking/GazeCal01_pupil.mat';

%% find sceneGeometry
sceneGeometryFileName = '~/Desktop/eyeTrackingDEMO/TOME_processing/session2_spatialStimuli/TOME_3020/050517/EyeTracking/GazeCal01_sceneGeometry.mat';
sceneGeometry = estimateSceneGeometry (pupilFileName,sceneGeometryFileName);

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
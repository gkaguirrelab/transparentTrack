% TEST_simulatedSceneGeometryConstraint

% This script generates a video composed of ellipses on the image plane
% that are pseudoPerspective projections of a model pupil. This video is
% then run through the pupil analysis pipeline, including derivation of the
% scene geometry and sceneConstrained ellipse fitting.

%% Clean up
clearvars
close all

% Reset the random number generator so that we have reproducible results
rng default

% 1. define synthetic data length
nFrames = 200;

% 2. define scene dimension in pixels
videoSizeX = 640;
videoSizeY = 480;

% Setup the video save 
sandboxDir = '~/Desktop/sceneGeometryFromPseudoPerspectiveTEST';

% check or make a directory for output
if exist(sandboxDir,'dir')==0
    mkdir(sandboxDir);
end

syntheticVideoName = fullfile(sandboxDir, 'synthetic.avi');

% create the video writer object
writerObj = VideoWriter(syntheticVideoName, 'Uncompressed AVI');
writerObj.FrameRate = 60;

% define positions in degrees of azimuth and elevation
allPupilAzi = randn(1,nFrames)*15; % in degrees
allPupilEle = randn(1,nFrames)*15; % in degrees

%% construct scene geometry
% we will construct a scene where the CoR sits right at the center, and the
% pupil rotates on a 150 pixel radius.

% define the eye sphere radius
eyeRadius =  125;
pupilRadius = 30;
sceneDistance = 1200;

% Set up the sceneGeometry
projectionModel = 'pseudoPerspective';
eyeCenter = [videoSizeX/2, videoSizeY/2, eyeRadius+sceneDistance];

emptyFrame = ones(videoSizeY, videoSizeX);
h = figure('visible','off');

% for each frame, create the projected ellipse and test the consistency of
% the returned theta values
for ii = 1:nFrames
    pupilAzi = allPupilAzi(ii); % in degrees
    pupilEle = allPupilEle(ii); % in degrees
    
    forwardProjectEllipseParams = pupilProjection_fwd(pupilAzi, pupilEle, pi*pupilRadius.^2, eyeCenter, eyeRadius, projectionModel);
    
    % make the plot and save it as a frame
    imshow(emptyFrame);
    hold on
    pFitImplicit = ellipse_ex2im(ellipse_transparent2ex(forwardProjectEllipseParams));
    fh=@(x,y) pFitImplicit(1).*x.^2 +pFitImplicit(2).*x.*y +pFitImplicit(3).*y.^2 +pFitImplicit(4).*x +pFitImplicit(5).*y +pFitImplicit(6);
    fimplicit(fh,[1, videoSizeX, 1, videoSizeY],'Color', 'b','LineWidth',1);
    axis equal
    axis off
    ylim([0 videoSizeY]);
    xlim([0 videoSizeX]);
    truesize;
    
    % get the frame
    tmpFrame = getframe(gca);
    [Yp, Xp] = ind2sub(size(squeeze(tmpFrame.cdata(:,:,1))),find(squeeze(tmpFrame.cdata(:,:,1))==0));
    
    % Fill the points
    fill(Xp,Yp,'k')
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
open(writerObj);

for ii=1:nFrames
    indexedFrame = rgb2ind(thisFrame(ii).cdata,cmap, 'nodither');
    
    % write video
    writeVideo(writerObj,indexedFrame);
end

close (writerObj);

%% use this video in the video pipeline 
perimeterFileName = fullfile(sandboxDir, 'synthetic_perimeter.mat');
pupilFileName = fullfile(sandboxDir, 'syntheticPerimeter_pupil.mat');
sceneGeometryFileName = fullfile(sandboxDir, 'syntheticPerimeter_sceneGeometry.mat');
sceneDiagnosticPlotFileName = fullfile(sandboxDir, 'syntheticPerimeter_sceneDiagnosticPlot.pdf');
finalFitVideoName = fullfile(sandboxDir, 'syntheticPerimeter_finalFit.avi');

findPupilPerimeter(syntheticVideoName,perimeterFileName,'verbosity','full','pupilRange', [20 400]);
fitPupilPerimeter(perimeterFileName, pupilFileName,'verbosity','full','ellipseTransparentLB',[],'ellipseTransparentUB',[],'nSplits',0);
sceneGeometry = estimateSceneGeometry(pupilFileName, sceneGeometryFileName,'sceneDiagnosticPlotFileName', sceneDiagnosticPlotFileName,'sceneDiagnosticPlotSizeXY', [videoSizeX videoSizeY], ...
    'projectionModel','pseudoPerspective','eyeRadius',eyeRadius, 'cameraDistanceInPixels',sceneDistance,'verbosity','full');
fitPupilPerimeter(perimeterFileName,pupilFileName,'sceneGeometryFileName',sceneGeometryFileName,'ellipseTransparentLB',[0, 0, 300, 0, 0],'ellipseTransparentUB',[videoSizeX,videoSizeY,20000,1.0, pi],'verbosity','full');
pupilData = smoothPupilArea(perimeterFileName, pupilFileName, sceneGeometryFileName,'verbosity','full');
makeFitVideo(syntheticVideoName, finalFitVideoName, 'pupilFileName',pupilFileName,'sceneGeometryFileName',sceneGeometryFileName,'perimeterFileName',perimeterFileName,'whichFieldToPlot','ellipseParamsAreaSmoothed_mean')


%% Verify that the scene geometry allows for the correct reconstruction of the eye position

ellipses = pupilData.ellipseParamsAreaSmoothed_mean;

for ii = 1:nFrames
    [reconstructedPupilAzi(ii), reconstructedPupilEle(ii), reconstructedPupilArea(ii)] = pupilProjection_inv(ellipses(ii,:),  [sceneGeometry.eyeCenter.X sceneGeometry.eyeCenter.Y, sceneGeometry.eyeCenter.Z], sceneGeometry.eyeRadius, sceneGeometry.meta.projectionModel);
end

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
plot(allPupilEle,reconstructedPupilEle, '.r')
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
rl = refline(0,pi*pupilRadius.^2);
rl.Color = 'k';
xlabel('frame')
ylabel('reconstructed pupil area')
axis square


%% Demonstrate how closely we have reconstructed the actual scene geometry
fprintf('Veridical scene geometry - eye center: [%0.1f, %0.1f, %0.1f], pupil radius: %0.1f \n',eyeCenter(1), eyeCenter(2), eyeCenter(3), eyeRadius);
fprintf('Estimated scene geometry - eye center: [%0.1f, %0.1f, %0.1f], pupil radius: %0.1f \n',sceneGeometry.eyeCenter.X, sceneGeometry.eyeCenter.Y, sceneGeometry.eyeCenter.Z, sceneGeometry.eyeRadius);


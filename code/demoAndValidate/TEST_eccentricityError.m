% TEST - full perspective eccentricity correction

% here, we use the blender eye model to test if the eccentricity correction
% formulas work as they should.

% first we generate a blender eye sweeping certain azimuth and elevations.
% We then track the eye using transparent track, in order to get the center
% position of the unconstrained ellipse fit (before even getting to the
% part where we reconstruct the scene geometry).
% We extract the ellipse center positions and apply a simplified version of
% Ahn's eccentricity correction formula and plot the results.


%% generate blender video

% Setup the video save location and path params

% This path should be defined in the eyeModelSupport local hook
codeDirectory = '~/Documents/MATLAB/toolboxes/eyemodelSupport/code';

% exportsDirectory = '~/Desktop/eccentricityCorrection60';
exportsDirectory = '~/Desktop/eccentricityCorrection';
pathParams.dataOutputDirFull = fullfile(exportsDirectory);
pathParams.dataSourceDirFull = fullfile(exportsDirectory);
pathParams.runName = 'simulation01';
videoName = fullfile(exportsDirectory, [pathParams.runName '_gray.avi']);

% check or make a directory for output
if exist(exportsDirectory,'dir')==0
    mkdir(exportsDirectory);
end


% features of the camera for the demo (derived from the blender model)
sensorSize = [36 24]; %mm
aperture = 2.0;
fieldOfViewDEG = 45;
sceneResolution = [640 480];
focalLengthPX = (sceneResolution(1)/2) / tand(45); % this is how it is defined in blender, but I am not sure it is right
focalLengthMM = 18;

focalLengthPX = (sceneResolution(2)/2) / tand(45); % this is the focal length that works

pixelSizeMM = sensorSize./sceneResolution;


% We find that objects in the image are 25% larger (in pixels) than what we
% would calculate here. For now, we apply this fudge factor while we try to
% understand how to calculate pixel size correctly.
fudgeFactor = 1.2;

pixelSizeMM = pixelSizeMM ./ fudgeFactor;


% eye radius (hardcoded in the python function)
eyeRadiusMM = 14;
eyeRadiusPX = eyeRadiusMM/pixelSizeMM(1);

% The camera distance is the distance of the camera from the center of
% rotation of the eye.
cameraDistanceMM = 50;

% The camera is hard-coded in the python routine to focus upon (i.e.,
% target) the front of the center of the eye. Therefore, the scene distance
% is equal to the camera distance + the eye raidus. We compute the scene
% distance in mm and then convert to pixels here.
sceneDistancePX = (cameraDistanceMM+eyeRadiusMM)/pixelSizeMM(1);


% Define the range of azimuth and elevation steps 
eleSteps = -20:20:0;
aziSweeps = -25:25:25;

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


% The eye is posed in the Blender model by defining the location in XYZ
% coordinate space of an eye target. The coordinate system that we use in
% the transparentTrack routines has the convention X = left/right; Y =
% down/up; Z = nearer / farther. The subsequent python routine transposes
% some of these dimensions to be appropriate for the Blender coordinate
% convention, but this need not concern us here.
%
% We calculate the target of the gaze using equations that describe the
% Fick axis rotation of the eye. In this system, the Y (down / up) position
% of the center of the pupil is indenpendent of the azimuthal rotation. The
% gazeTargetPositions that we generate correspond to a point that is
% directly on the surface of the eye in the center of the pupil.
gazeTargetPositionX = (eyeRadiusMM).*sind(allPupilAzi).*cosd(allPupilEle);
gazeTargetPositionY = (eyeRadiusMM).*sind(allPupilEle);
gazeTargetPositionZ = -(eyeRadiusMM).*cosd(allPupilAzi).*cosd(allPupilEle);

% Define the pupil radius in pixels and the degree of eye closedness
pupilRadiusMM = 2*ones(1,length(gazeTargetPositionX));
pupilRadiusPX = pupilRadiusMM/pixelSizeMM(1);
eyeClosedness = 0*ones(1,length(gazeTargetPositionX));


%% generate eye movie
tbUse('eyemodelSupport', 'reset','full')

generateEyeMovie(codeDirectory,exportsDirectory,gazeTargetPositionX, gazeTargetPositionY, gazeTargetPositionZ, pupilRadiusMM, eyeClosedness,cameraDistanceMM)

% rename file
movefile(fullfile(exportsDirectory,'pupil_movie.avi'),fullfile(exportsDirectory,[pathParams.runName '_gray.avi']));

tbUse('transparentTrack', 'reset','full')
%% run tracking pipeline up to ellipse unconstrained fit.
runVideoPipeline( pathParams, ...
    'verbosity', 'full', 'useParallel',false, 'catchErrors', false,...
    'pupilFrameMask', [60 60], 'maskBox', [0.9 0.9], 'pupilGammaCorrection',0.7, 'pupilRange', [30 80], ...
    'overwriteControlFile',true, 'glintPatchRadius', 10, ...
    'ellipseTransparentLB',[0,0, 20, 0, 0],...
    'ellipseTransparentUB',[sceneResolution(1),sceneResolution(2), 20000, 1.0, pi],...
    'candidateThetas',0:pi/16:2*pi,...    
    'badFrameErrorThreshold', 8, ...
    'sceneGeometryLB',[0, 0, sceneDistancePX, 25],'sceneGeometryUB',[640, 480, sceneDistancePX, 500],...
    'skipStageByNumber',1, ...
    'lastStageByNumber', 6, 'makeFitVideoByNumber', 6);

%% compare ellipse centers to the position in which we would expect to see the pupil centers.

% get the ellipse data
load(fullfile(pathParams.dataOutputDirFull,[pathParams.runName '_pupil.mat']));

% this is where the ellipse center is tracked
ellipseXpos = pupilData.ellipseParamsUnconstrained_mean(:,1);
ellipseYpos = pupilData.ellipseParamsUnconstrained_mean(:,2);



% ellipse center, circle center and error formulas
for ii = 1:length(allPupilAzi)

    centerOfCircleX(ii) = (sceneResolution(1)/2) + eyeRadiusPX*sind(allPupilAzi(ii))*cosd(-allPupilEle(ii));
    centerOfCircleY(ii) = (sceneResolution(2)/2) + eyeRadiusPX*sind(-allPupilEle(ii));
    centerOfCircleZ(ii) = eyeRadiusPX - eyeRadiusPX*cosd(allPupilAzi(ii))*cosd(-allPupilEle(ii));
    
    pupilCenterToCameraDistancePx(ii) = sceneDistancePX+centerOfCircleZ(ii);

% 
% centerOfCircleX(ii) = (sceneResolution(1)/2) + (focalLengthPX *tand(allPupilAzi(ii)));
% 
centerOfEllipseX(ii) = (sceneResolution(1)/2) + ((focalLengthPX *((pupilCenterToCameraDistancePx(ii)^2 * sind(allPupilAzi(ii)) * cosd(allPupilAzi(ii))) + (pupilRadiusPX(1)^2 * sind(allPupilAzi(ii)) * cosd(allPupilAzi(ii)))) / ...
     ((pupilCenterToCameraDistancePx(ii)^2 * cosd(allPupilAzi(ii))) - (pupilRadiusPX(1)^2 *(sind(allPupilAzi(ii)))^2)))*fudgeFactor);
 
 
 
 translationComponentX(ii) = ( (focalLengthPX * (eyeRadiusPX*sind(-allPupilEle(ii))/pupilCenterToCameraDistancePx(ii)) * (sind(allPupilAzi(ii)))^2) / ...
     ((pupilCenterToCameraDistancePx(ii)/pupilRadiusPX(ii))^2 - (sind(allPupilAzi(ii)))^2 ) ).*fudgeFactor;
 
 geometricCenterErrorX(ii) = ( (focalLengthPX * sind(allPupilAzi(ii))) / ...
     (  (pupilCenterToCameraDistancePx(ii)/pupilRadiusPX(ii))^2 - (sind(allPupilAzi(ii)))^2 ) ).*fudgeFactor;
 
 
 

 % 
% centerOfCircleY(ii) = (sceneResolution(2)/2) + (-focalLengthPX *tand(allPupilEle(ii)));
% 
 centerOfEllipseY(ii) = (sceneResolution(2)/2) + ((-focalLengthPX *((pupilCenterToCameraDistancePx(ii)^2 * sind(allPupilEle(ii)) * cosd(allPupilEle(ii))) + (pupilRadiusPX(1)^2 * sind(allPupilEle(ii)) * cosd(allPupilEle(ii)))) / ...
     ((pupilCenterToCameraDistancePx(ii)^2 * cosd(allPupilEle(ii))) - (pupilRadiusPX(1)^2 *(sind(allPupilEle(ii)))^2)))*fudgeFactor);
% 

 translationComponentY(ii) = ( (focalLengthPX * (eyeRadiusPX*sind(allPupilAzi(ii))*cosd(allPupilEle(ii))/pupilCenterToCameraDistancePx(ii)) * (sind(allPupilEle(ii)))^2) / ...
     ((pupilCenterToCameraDistancePx(ii)/pupilRadiusPX(ii))^2 - (sind(allPupilEle(ii)))^2 ) ).*fudgeFactor;
 
 geometricCenterErrorY(ii) = ((focalLengthPX*sind(allPupilEle(ii)))/ ...
     ((pupilCenterToCameraDistancePx(ii)/pupilRadiusPX(1))^2 - (sind(allPupilEle(ii)))^2)).*fudgeFactor;
end

% compare centers of ellipses (the tracked ones should more or less overlay
% where we estimate the ellipses centers to be).
% figure
% subplot(2,1,1)
% plot(centerOfCircleX)
% hold on
% plot(centerOfEllipseX)
% hold on
% plot (ellipseXpos)
% 
% xlabel('Frames')
% ylabel('Horizontal position in pixels')
% legend ('theoretical center of circle', 'anh','center of fitted ellipse that we measure')
% title ('Teoretical and tracked X position of ellipse center')
% 
% subplot(2,1,2)
% plot(centerOfCircleY)
% hold on
% plot(centerOfEllipseY)
% hold on
% plot (ellipseYpos+1)
% xlabel('Frames')
% ylabel('Vertical position in pixels')
% legend ('theoretical center of circle',  'anh','center of fitted ellipse that we measure')
% title ('Teoretical and tracked Y position of ellipse center')

figure
plot(centerOfCircleX,centerOfCircleY)
hold on
plot(centerOfEllipseX,centerOfEllipseY)
plot (ellipseXpos,ellipseYpos+1)
plot (centerOfCircleX+geometricCenterErrorX,centerOfCircleY+geometricCenterErrorY)

xlim ([0 sceneResolution(1)])
ylim ([0 sceneResolution(2)])
set(gca,'Ydir','reverse')
legend ('our forward model',  'anh','center of fitted ellipse that we measure','our forward model with Anh geometric center fix')
title('Distortion correction of pupil center pattern on screen')

% %% get the perspective projection fwd to put the ellipses centers where we expect them to be
% pupilArea = pi .* pupilRadiusPX.^2;
% eyeCenter = [320 240 eyeRadiusPX+sceneDistancePX];
% 
% for ii = 1:length(allPupilAzi)
%     reconstructedTransparentEllipse(ii,:) = pupilProjection_fwd(allPupilAzi(ii), allPupilEle(ii), pupilArea(ii), eyeCenter, eyeRadiusPX, 'pseudoPerspective');
% end
% 
% projectedCentersX = reconstructedTransparentEllipse(:,1);
% projectedCentersY = reconstructedTransparentEllipse(:,2);
% 
% figure
% plot (ellipseXpos,ellipseYpos+1)
% hold on
% plot(projectedCentersX,projectedCentersY)
% legend ('blender', 'fwd projection')
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
exportsDirectory = '~/Desktop/blenderEyeModelSimulation';
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


% eye radius and camera distance (hardcoded in the python function)
eyeRadiusMM = 14;
cameraDistanceMM = 60;

eyeRadiusPX = eyeRadiusMM/pixelSizeMM(1);
sceneDistancePX = cameraDistanceMM/pixelSizeMM(1);


%define azimuth and elevation timeseries
% Define the range of azimuth and elevation steps 
eleSteps = -20:5:20;
aziSweeps = -30:3:30;

% shorter versions here. Decomment the ones needed
% ORIZONTAL SWEEP ONLY
% eleSteps = 0;
% aziSweeps = -30: 5:30;

% VERTICAL STEP ONLY
% eleSteps = -10: 2: 10;
% aziSweeps = 0;

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

rotationArm = abs(eyeRadiusMM);
% pupil position based on the angles (eye target)
pupilXpos = (rotationArm).*sind(allPupilAzi).*cosd(allPupilEle);
pupilYpos = (rotationArm).*sind(allPupilEle);
pupilZpos =-(rotationArm).*cosd(allPupilAzi).*cosd(allPupilEle);

pupilRadiusMM = 2*ones(1,length(pupilXpos));
pupilRadiusPX = pupilRadiusMM/pixelSizeMM(1);
eyeClosedness = 0*ones(1,length(pupilXpos));


%% generate eye movie
tbUse('eyemodelSupport', 'reset','full')

generateEyeMovie(codeDirectory,exportsDirectory,pupilXpos, pupilYpos, pupilZpos, pupilRadiusMM, eyeClosedness,cameraDistanceMM)

% rename file
movefile(fullfile(exportsDirectory,'pupil_movie.avi'),fullfile(exportsDirectory,[pathParams.runName '_gray.avi']));

tbUse('transparentTrack', 'reset','full')
%% run tracking pipeline up to ellipse unconstrained fit.
runVideoPipeline( pathParams, ...
    'verbosity', 'full', 'useParallel',false, 'catchErrors', false,...
    'pupilFrameMask', [60 60], 'maskBox', [0.9 0.9], 'pupilGammaCorrection',0.7, 'pupilRange', [10 60], ...
    'overwriteControlFile',true, 'glintPatchRadius', 10, ...
    'projectionModel', 'orthogonal', ...
    'ellipseTransparentLB',[0,0, 20, 0, 0],...
    'ellipseTransparentUB',[sceneResolution(1),sceneResolution(2), 20000, 1.0, pi],...
    'candidateThetas',0:pi/16:2*pi,...
    'badFrameErrorThreshold', 8, ...
    'skipStageByNumber',1, ...
    'lastStageByNumber', 6, 'makeFitVideoByNumber', 6);

%% compare ellipse centers to the position in which we would expect to see the pupil centers.

% get the ellipse data
load(fullfile(pathParams.dataOutputDirFull,[pathParams.runName '_pupil.mat']));

% this is where the ellipse center is tracked
ellipseXpos = pupilData.ellipseParamsUnconstrained_mean(:,1);
ellipseYpos = pupilData.ellipseParamsUnconstrained_mean(:,2);


planeTiltAngleX = allPupilAzi;
planeTiltAngleY = allPupilEle;


l = (sceneDistancePX+focalLengthPX) .*(cosd(abs(planeTiltAngleY)));
cameraDistancePX  = (sceneDistancePX+focalLengthPX);

% ellipse center, circle center and error formulas
for ii = 1:length(planeTiltAngleX)

centerOfEllipseX(ii) = (sceneResolution(1)/2) + (focalLengthPX *((sceneDistancePX^2 * sind(planeTiltAngleX(ii)) * cosd(planeTiltAngleX(ii))) + (pupilRadiusPX(1)^2 * sind(planeTiltAngleX(ii)) * cosd(planeTiltAngleX(ii)))) / ...
    ((sceneDistancePX^2 * cosd(planeTiltAngleX(ii))) - (pupilRadiusPX(1)^2 *(sind(planeTiltAngleX(ii)))^2)));

reconstructedErrorX(ii) = (focalLengthPX*((cameraDistancePX)/l(ii))*sind(planeTiltAngleX(ii))*cosd(planeTiltAngleX(ii)))/ ...
    ((l(ii)/pupilRadiusPX(1))^2 - (sind(planeTiltAngleX(ii)))^2);    

centerOfEllipseY(ii) = (sceneResolution(2)/2) + (-focalLengthPX *((sceneDistancePX^2 * sind(planeTiltAngleY(ii)) * cosd(planeTiltAngleY(ii))) + (pupilRadiusPX(1)^2 * sind(planeTiltAngleY(ii)) * cosd(planeTiltAngleY(ii)))) / ...
    ((sceneDistancePX^2 * cosd(planeTiltAngleY(ii))) - (pupilRadiusPX(1)^2 *(sind(planeTiltAngleY(ii)))^2)));

reconstructedErrorY(ii) = (focalLengthPX*((cameraDistancePX)/l(ii))*sind(planeTiltAngleY(ii))*cosd(planeTiltAngleY(ii)))/ ...
    ((l(ii)/pupilRadiusPX(1))^2 - (sind(planeTiltAngleY(ii)))^2);
end

% compare centers of ellipses (the tracked ones should more or less overlay
% where we estimate the ellipses centers to be).
figure
subplot(2,1,1)
plot(centerOfEllipseX)
hold on
plot (ellipseXpos)
hold on
plot (ellipseXpos - reconstructedErrorX')
xlabel('Frames')
ylabel('Horizontal position in pixels')
legend ('theoretical center of ellipse according to Ahn formula', 'center of fitted ellipse that we measure', 'theoretical ellipse center - eccentricity error (Ahn)')
title ('Teoretical and tracked X position of ellipse center')

subplot(2,1,2)
plot(centerOfEllipseY)
hold on
plot (ellipseYpos + 1) % adding 1 for indexing discrepancy
hold on
plot ((ellipseYpos + 1) - reconstructedErrorY')
xlabel('Frames')
ylabel('Vertical position in pixels')
legend ('theoretical center of ellipse according to Ahn formula', 'center of fitted ellipse that we measure', 'theoretical ellipse center - eccentricity error (Ahn)')
title ('Teoretical and tracked Y position of ellipse center')

figure
plot(centerOfEllipseX,centerOfEllipseY)
hold on
plot (ellipseXpos,ellipseYpos+1)
hold on
plot (ellipseXpos - reconstructedErrorX',(ellipseYpos + 1) - reconstructedErrorY')
xlim ([0 sceneResolution(1)])
ylim ([0 sceneResolution(2)])
set(gca,'Ydir','reverse')
legend ('theoretical center of ellipse according to Ahn formula', 'center of fitted ellipse that we measure', 'theoretical ellipse center - eccentricity error (Ahn)')
title('Distortion correction of pupil center pattern on screen')

% 
% % plot the measured errors
% figure
% subplot(2,1,1)
% plot (measuredErrorX)
% ylabel ('Error in PX')
% xlabel ('Frames')
% title('Ellipse centerX - projected circle centerX')
% subplot(2,1,2)
% plot (measuredErrorY)
% ylabel ('Error in PX')
% xlabel ('Frames')
% title('Ellipse centerY - projected circle centerY')
% % note that in the pupilData there is a big fitting error on the first ellipse in the Y
% % dimension




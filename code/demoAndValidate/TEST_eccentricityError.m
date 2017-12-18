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

exportsDirectory = '~/Desktop/eccentricityCorrection60';
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
% focalLengthPX = (focalLengthMM/sensorSize(1))*sceneResolution(1); % this is the formula I would use

pixelSizeMM = sensorSize./sceneResolution;


% eye radius and camera distance (hardcoded in the python function)
eyeRadiusMM = 14;
cameraDistanceMM = 60;

eyeRadiusPX = eyeRadiusMM/pixelSizeMM(1);
sceneDistancePX = cameraDistanceMM/pixelSizeMM(1);


%define azimuth and elevation timeseries
% Define the range of azimuth and elevation steps 
% eleSteps = -20:5:20;
% aziSweeps = -30:3:30;

% shorter version
eleSteps = [0];
aziSweeps = [-30: 5:30];

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


planeTiltAngle = -30:5:30; % since we only do an Azi sweep.

planeTiltAngle = linspace(-22, 22, 13); % looks like this apparent angle works. WHY?

l = (sceneDistancePX+focalLengthPX) .*(cosd(abs(planeTiltAngle)));
cameraDistancePX  = (sceneDistancePX+focalLengthPX);
% ellipse center position forumlas
for ii = 1:size(planeTiltAngle,2)

centerOfEllipse(ii) = focalLengthPX *((sceneDistancePX^2 * sind(planeTiltAngle(ii)) * cosd(planeTiltAngle(ii))) + (pupilRadiusPX(1)^2 * sind(planeTiltAngle(ii)) * cosd(planeTiltAngle(ii)))) / ...
    ((sceneDistancePX^2 * cosd(planeTiltAngle(ii))) - (pupilRadiusPX(1)^2 *(sind(planeTiltAngle(ii)))^2)) + 320;

reconstructedErrorX(ii) = (focalLengthPX*((cameraDistancePX)/l(ii))*sind(planeTiltAngle(ii))*cosd(planeTiltAngle(ii)))/ ...
    ((l(ii)/pupilRadiusPX(1))^2 - (sind(planeTiltAngle(ii)))^2)
end

% compare centers of ellipses (the tracked ones should more or less overlay
% where we estimate the ellipses centers to be).
figure
plot(centerOfEllipse)
hold on
plot (ellipseXpos)
legend ('center of ellipse according to Ahn formula', 'center of fitted ellipse that we measure')



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



%% %% THIS PART NEEDS MORE THINKING :(
% 


% % pupil center position on scene (using pinhole model)
% % this is where the pupil center would be with no distortion.
% pupilXScene = (sceneResolution(1)/2) + ((focalLengthPX/sceneDistancePX) .*(pupilXpos'./pixelSizeMM(1))) ;
% pupilYScene = (sceneResolution(2)/2) + ((focalLengthPX/sceneDistancePX) .*(pupilYpos'./pixelSizeMM(1))) ;




% % Ahn's estimation of the error (OBJECT PLANE)
% planeTiltAngle = allPupilAzi; % since we only do an Azi sweep.
% 
% % 1. project back the center
% ellipseXObject = (ellipseXpos - sceneResolution(1)/2) .* sceneDistancePX/focalLengthPX;
% ellipseYObject = -focalLengthPX;
% 
% %2. calculate error from data
% errorOnObject = [ellipseXObject - pupil3DPX(:,1) ellipseYObject - pupil3DPX(:,2)] ;
% 
% % calculate error from Ahn's formula
% AhnErrorY = 0;
% AhnErrorX = (-pupilRadiusPX.^2 .* sind (allPupilAzi))./sceneDistancePX;
% 
% 
% 
% % get tilt angles from PX values
% centralPlaneN = [0 0 eyeRadiusPX];
% pupilCenter3D = -[pupilXpos' pupilYpos' pupilZpos']; %minus because of blender convention
% for ii = 1: size(pupilCenter3D, 1)
%     planeTiltAngle(ii) = atan2(norm(cross(centralPlaneN,pupilCenter3D(ii,:))),dot(centralPlaneN,pupilCenter3D(ii,:)));
% end
%  
% 
% % get the errors
% alfa = deg2rad(fieldOfViewDEG); % field of view
% for ii = 1: size(pupilCenter3D, 1)
% l(ii) = sceneDistancePX*(cos(planeTiltAngle(ii) - alfa));
% reconstructedErrorX(ii) = (focalLengthPX*(pupil3DPX(ii,1)/l(ii))*(sin(planeTiltAngle(ii)))^2) / ((l(ii)/pupilRadiusPX(ii))^2 - (sin(planeTiltAngle(ii)))^2);
% reconstructedErrorY(ii) = (-focalLengthPX*(sceneDistancePX/l(ii))*sin(planeTiltAngle(ii))*cos(planeTiltAngle(ii)))/ ((l(ii)/pupilRadiusPX(ii))^2 - (sin(planeTiltAngle(ii)))^2);
% 
% reconstructedError(ii) = (-pupilRadiusPX(ii)^2 * sin(planeTiltAngle(ii))) /l(ii);
% end
% %% plot veridical center locations and corrected ellipses center locations
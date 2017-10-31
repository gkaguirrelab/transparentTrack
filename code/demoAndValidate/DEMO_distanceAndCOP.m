% DEMO find distances and COP

%% housekeeping
close all
clear all
%clc

%% load pupil data
pupilFileName = '~/Desktop/eyeTrackingDEMO/TOME_processing/session2_spatialStimuli/TOME_3020/050517/EyeTracking/GazeCal01_pupil.mat';
sizeCalibrationFileName = '~/Desktop/eyeResponseDemo/TOME_processing/session2_spatialStimuli/TOME_3020/050517/EyeTracking/tfMRI_RETINO_PA_run01_sizeCalFactors.mat';

load(pupilFileName)
ellipses = pupilData.pPosteriorMeanTransparent(1:500,:);

load(sizeCalibrationFileName)
averageEyeballMM = 50; % mm
eyeballRadiusInImagePlanePixels  = averageEyeballMM * sizeCalFactors.horizontalPxPerMm * 2;

%% find the most circular ellipse and use the X Y coordinate as the sceneCOP
[minEccentricity, minEccentricityIDX] = min(pupilData.pPosteriorMeanTransparent(:,4));

% In the orthogonal projection case, the distance of the center of eyeball
% from the image plane is equivalent to the radius of the eyeball

%% first guess for the CoP is the location of the most circular ellipse
x0 = [pupilData.pPosteriorMeanTransparent(minEccentricityIDX,1) ...
    pupilData.pPosteriorMeanTransparent(minEccentricityIDX,2) ...
    eyeballRadiusInImagePlanePixels];

errorWeights = (1./pupilData.fitError)';

% define an anonymous function to measure SSQ error 
errorFunc = @(x) sqrt(nansum(errorWeights.*distanceToCandidateEyeBallCenter(ellipses,x).^2));
bestFitCoP = fmincon(errorFunc, x0);

[distances, ellipsesCOPs] = distanceToCandidateEyeBallCenter(ellipses,bestFitCoP);

%% plot all the COPs

figure
plot(ellipsesCOPs(:,1), ellipsesCOPs(:,2), '.')
hold on
plot(x0(1),x0(2), 'xr')
plot(bestFitCoP(1),bestFitCoP(2), 'og')
title('Estimate Position of the Center of Projeciton for fitted ellipses')
legend('Ellipses COP', 'COP of the most circular ellipse (reference COP)','Best fit CoP')
xlim([0 320] * 2)
ylim([0 240] * 2)
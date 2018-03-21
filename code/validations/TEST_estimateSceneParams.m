%% TEST_estimateSceneParams
% Examine the ability of the routines to estimate an unknown scene geometry
%
% Description:
%   A core operation of transparentTrack is estimation of the extrinsic
%   camera translation vector from a set of pupil ellipses on the image
%   plane. The accuracy of this estimation depends in part upon the
%   accuracy of specification of the parameters of eye, and-most
%   critically--upon the axial length of the eye and its rotation center.
%   Here, we examine the ability of the scene estimation routine to recover
%   the translation vector given a set of ellipses generated from a
%   veridical model.
%
%   We first observe that the routine finds the X and Y position accurately
%   regardless of Z position. We then examine the ability to recover depth.
%   We perform the calculation assuming the default eye model parameters,
%   but generate the ellipses using modeled eyes with axial lengths that
%   deviate from the assumed value.
%
%   The axial length of the eye varies considerably across individuals. The
%   precise length may be measured (e.g., using the Zeiss IOL Master). If
%   such a measurement is available, it may be used to inform the eye model
%   and produce a more accurate estimation of scene geometry. If this is
%   unavailable, the investigator may have available the spherical
%   refractive error of the subject (i.e., diopters of correction for
%   myopia or hyperopia). Axial length is correlated with spectacle
%   correction, but not perfectly so. In our simulation, we consider the
%   case in which we have the spectacle correction for the subject and ask
%   how large an error in estimation of the axial length may be expected.
%   To compute this, we examined the data from the Reykjavik Eye Study:
%
%       Olsen, Thomas, et al. "On the ocular refractive components: the
%       Reykjavik Eye Study." Acta Ophthalmologica 85.4 (2007): 361-366.
%
%   Figures 1 and 4 give the observed distribution of refractive error and
%   axial lengths in their population of 723 adult, right eyes. I digitized
%   the values in the plots, fit the data with Gaussian distribution and
%   retained the sigma values. The paper reports that, across eyes, there
%   was a correlation of -0.59 between axial length and spherical
%   refractive error. Given this information, we can calculate the width of
%   the distribution of axial lengths given the refractive error of a
%   subject. This calculation yields a sigma of 0.9384 mm. In our examination
%   of the sceneGeometry estimation, therefore, we examine the error in the
%   estimation of the extrinsic translation vector given an error in
%   specification of axial length of 1SD (+- 0.9384 mm) or 2SDs.
%


% Set up a destination to save the files as they are generated one by one.
% We save the data as they are generated because this can take a while and
% we'd like to recover gracefully from an interruption of the process.
thisComputer = computer;
switch thisComputer
    case 'MACI64'
        outputFileStem = fullfile('~','Dropbox (Aguirre-Brainard Lab)','TOME_analysis','gka_simulationTests','TEST_estimateSceneParams','sceneGeometry');
    case 'GLNXA64'
        outputFileStem = fullfile('~','TEST_estimateSceneParams','sceneGeometry');
end


%% Analyze data extracted from Olsen 2007
% Histogram data from Figure 1 of Olsen 2007
res_refractionVals = [-12	-7	-6	-5	-4	-3	-2	-1	0	1	2	3	4	5	6	7	8];
res_refractionCounts = [1	1	4	5	25	19	28	79	149	189	114	60	23	6	7	5	1];
tmpFit = fit(res_refractionVals',res_refractionCounts','gauss1');
sigmaRefractionDiopters  = tmpFit.c1;

% Histogram data from Figure 4 of Olsen 2007
res_axialLengthVals = [21 22 23 24 25 26 27 28];
res_axialLengthCounts = [7	40	199	291	126	51	5	2];
tmpFit = fit(res_axialLengthVals',res_axialLengthCounts','gauss1');
sigmaLengthMm = tmpFit.c1;

% Reported correlation between axial length and diopters of refractive
% error in Olsen 2007
p = -0.59;

% Calculate the width of the distribution of axial lengths given the
% spherical refraction. This draws upon the conditional expectation of a
% X given Y in a bivariate normal distribution
conditionalSigmaLength = sqrt((1-p^2)*sigmaLengthMm);

% Obtain the axial length of the default model eye.
defaultSceneGeometry = createSceneGeometry();
defaultAxialLength = defaultSceneGeometry.eye.axialLength;

% Create and save the default scene geometry
outputFile = [outputFileStem '_default.mat'];
save(outputFile,'defaultSceneGeometry');

%% Create a veridical sceneGeometry with the specified axial lenght
axialLength = defaultAxialLength + (axialErrorMultiplier * conditionalSigmaLength);

veridicalSceneGeometry = createSceneGeometry('axialLength',axialLength);
outputFile = [outputFileStem '_veridical_axialLength=' num2str(axialLength,'%2.2f') '.mat'];
save(outputFile,'veridicalSceneGeometry');

% Assemble the ray tracing functions
rayTraceFuncs = assembleRayTraceFuncs( veridicalSceneGeometry );

% Create a set of ellipses using the veridical geometry and randomly
% varying pupil radii.
rng('default')
ellipseIdx=1;
for azi=-15:15:15
    for ele=-15:15:15
        eyePose=[azi, ele, 0, 2+(randn()./5)];
        pupilData.initial.ellipses.values(ellipseIdx,:) = pupilProjection_fwd(eyePose, veridicalSceneGeometry, rayTraceFuncs);
        pupilData.initial.ellipses.RMSE(ellipseIdx,:) = 1;
        ellipseIdx=ellipseIdx+1;
    end
end

startTime=datetime('now');
result = estimateSceneParams(pupilData,'','axialLength',defaultAxialLength, ...
    'useParallel',false,'verbosity','full','ellipseArrayList',1:1:ellipseIdx-1, ....
    'nBADSsearches',100,'useRayTracing',false, ...
	'sceneParamsLB',[-20; -20; 40; 1; 1], ...
	'sceneParamsUB',[20; 20; 200; 1; 1], ...
	'sceneParamsLBp',[-5; -5; 80; 1; 1], ...
	'sceneParamsUBp',[5; 5; 160; 1; 1]);

endTime=datetime('now');
result.startTime = startTime;
result.endTime = endTime;
outputFile = [outputFileStem '_result_axialLength=' num2str(axialLength,'%2.2f') '.mat'];
save(outputFile,'result');


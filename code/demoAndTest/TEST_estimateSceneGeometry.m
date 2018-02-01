%% TEST_estimateSceneGeometry
% Examine the ability of the routines to estimate an unknown scene geometry
%
% Description:
%   A core operation of transparentTrack is estimation of the extrinsic
%   camera translation vector from a set of pupil ellipses on the image
%   plane. The accuracy of this estimation depends in part upon the
%   accuracy of specification of the parameters of eye, and?-most
%   critically--upon the axial length of the eye and its rotation center.
%   Here, we examine the ability of the scene estimation routine to recover
%   the translation vector given a set of ellipses generated from a
%   veridical model. The search is conducted for a range of X, Y, and Z
%   translation values to demonstrate the generality of the accuracy of the
%   solution. Additionally, we perform the calculation assuming the default
%   eye model parameters, but generate the ellipses using modeled eyes with
%   axial lengths that deviate from the assumed value.
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
%   axial lengths in their population of 723 adult, right eyes. I digitzed
%   the values in the plots, fit the data with Gaussian distribution amnd
%   retained the sigma values. The paper reports that, across eyes, there
%   was a correlation of -0.59 between axial length and spherical
%   refractive error. Given this information, we can calculate the width of
%   the distribution of axial lengths given the refractive error of a
%   subject. This calculation yields a sigma of 0.9384 mm. In our examination
%   of the sceneGeometry estimation, therefore, we examine the error in the
%   estimation of the extrinsic translation vector given an error in
%   specification of axial length of 1SD (+- 0.9384 mm) or 2SDs.
%

% Clean up for the simulation
clear all
close all

% Set up a destination to save the files as they are generated one by one.
% We save the data as they are generated because this can take a while and
% we'd like to recover gracefully from an interuption of the process.
thisComputer = computer;
switch thisComputer
    case 'MACI64'
        outputFileStem = fullfile('~','Dropbox (Aguirre-Brainard Lab)','TOME_analysis','gka_simulationTests','TEST_estimateSceneGeometry');
    case 'GLNXA64'
        outputFileStem = fullfile('~','TEST_estimateSceneGeometry');
end

%% Analyze data extracted from Olsen 2007
% Histogram data from Figure 1 of Olsen 2007
res_refractionVals = [-12	-7	-6	-5	-4	-3	-2	-1	0	1	2	3	4	5	6	7	8];
res_refractionCounts = [1	1	4	5	25	19	28	79	149	189	114	60	23	6	7	5	1];
tmpFit = fit(res_refractionVals',(res_refractionCounts./max(res_refractionCounts))','gauss1');
sigmaRefractionDiopters  = tmpFit.c1;

% Histogram data from Figure 4 of Olsen 2007
res_axialLengthVals = [21 22 23 24 25 26 27 28];
res_axialLengthCounts = [7	40	199	291	126	51	5	2];
tmpFit = fit(res_axialLengthVals',(res_axialLengthCounts./max(res_axialLengthCounts))','gauss1');
sigmaLengthMm = tmpFit.c1;

% Reported correlation between axial length and diopters of refractive
% error in Olsen 2007
p = -0.59;

% Calculate the width of the distribution of axial lengths given the
% spherical refraction. This draws upon the conditional expectation of a
% X given Y in a bivariate normal distribution
conditionalSigmaLength = sqrt((1-p^2)*sigmaLengthMm);

% Obtain the axial length of the default model eye. The axial length is
% given by the depth of the center of the posterior chamber, plus the
% radius of the axial dimension of the posterior chamber
defaultSceneGeometry = estimateSceneGeometry([],[]);
defaultAxialLength = -(defaultSceneGeometry.eye.posteriorChamberCenter(1)-defaultSceneGeometry.eye.posteriorChamberRadii(1));


%% Run the simulation
% Create a set of 9 pupil radii which will be used for all the the
% searches. We use random values to demonstrate that the solution is not
% dependent upon the area of the pupil ellipse, just position and shape
pupilRadii = 2+(randn(9,1)./5);

% Loop over simulation conditions and estimate scene geometry
resultIdx = 1;
for axialErrorMultiplier = -2:1:2
    for cameraX = -5:5:5
        for cameraY = -5:5:5
            for cameraZ = 125:25:175
                % We will create the ellipses using an axial length for the
                % eye that can is off by + and - 2SD of the conditional
                % distribution of axial lengths given knowledge of the
                % spherical refraction of the eye
                axialLength = defaultAxialLength + (axialErrorMultiplier * conditionalSigmaLength);
                veridicalSceneGeometry = estimateSceneGeometry([],[],'eyeLaterality','Right','axialLength',axialLength);
                veridicalSceneGeometry.extrinsicTranslationVector = [cameraX; cameraY; cameraZ];
                
                % Assemble the ray tracing functions
                rayTraceFuncs = assembleRayTraceFuncs( veridicalSceneGeometry );
                
                % Create a set of ellipses from the veridial geometry
                ellipseIdx=1;
                for azi=-15:15:15
                    for ele=-15:15:15
                        eyePoses=[azi, ele, 0, pupilRadii(ellipseIdx)];
                        pupilData.initial.ellipses.values(ellipseIdx,:) = pupilProjection_fwd(eyePoses, veridicalSceneGeometry, rayTraceFuncs);
                        pupilData.initial.ellipses.RMSE(ellipseIdx,:) = 1;
                        ellipseIdx=ellipseIdx+1;
                    end
                end
                
                % Estimate the scene Geometry
                estimatedSceneGeometry = estimateSceneGeometry(pupilData,'','useParallel',true,'nWorkers',8,'ellipseArrayList',1:1:ellipseIdx-1);

                % Save the veridical and estimated results
                outputFile = [outputFileStem '_vsg_' num2str(resultIdx) '.mat'];
                save(outputFile,'veridicalSceneGeometry');
                outputFile = [outputFileStem '_esg_' num2str(resultIdx) '.mat'];
                save(outputFile,'estimatedSceneGeometry');

                % Iterate the result index
                resultIdx = resultIdx+1;
            end
        end
    end
end


%% Reload the data and plot
figure
plotColors=[1 0 0; 0.5 0 0; 0 0 0; 0 0.5 0; 0 1 0];
resultIdx = 1;
for axialErrorMultiplier = -2:1:2
    for cameraX = -5:5:5
        for cameraY = -5:5:5
            for cameraZ = 125:25:175

                inputFile = [outputFileStem '_vsg_' num2str(resultIdx) '.mat'];
                load(inputFile);

                plot3([cameraX vsg.extrinsicTranslationVector(1)],...
                    [cameraY vsg.extrinsicTranslationVector(2)],...
                    [cameraZ vsg.extrinsicTranslationVector(3)],...
                    'Color',plotColors(axialErrorMultiplier+3,:),...
                    'Marker','o');
                
                hold on
            end
        end
    end
end                

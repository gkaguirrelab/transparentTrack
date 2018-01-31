%% TEST_estimateSceneGeometry
% Examine the ability of the routines to estimate an unknown scene geometry
%
% Description:

%
clear all
close all

outputFileStem = fullfile('~','Dropbox (Aguirre-Brainard Lab)','TOME_analysis','gka_simulationTests','TEST_estimateSceneGeometry');

% The Reykjavik Eye Study indicates that axial length is related to
% spherical refractive error with a correlation of -0.59, and a slope of
% -1.24, meaning that there was a decrease of 1.24 in diopters of spherical
% refractive error for every mm lengthening of the eye in their population
% of 723 right eyes. Overall, knowing the
%
%   Olsen, Thomas, et al. "On the ocular refractive components: the
%   Reykjavik Eye Study." Acta Ophthalmologica 85.4 (2007): 361-366.
%
% to determine what the 2 SD range of error in prediction of axial length
% is given only spherical refraction of the studied subject

res_refractionVals = [-12	-7	-6	-5	-4	-3	-2	-1	0	1	2	3	4	5	6	7	8];
res_refractionCounts = [1	1	4	5	25	19	28	79	149	189	114	60	23	6	7	5	1];
refractionGauss = fit(res_refractionVals',(res_refractionCounts./max(res_refractionCounts))','gauss1');

res_axialLengthVals = [21 22 23 24 25 26 27 28];
res_axialLengthCounts = [7	40	199	291	126	51	5	2];
lengthGauss = fit(res_axialLengthVals',(res_axialLengthCounts./max(res_axialLengthCounts))','gauss1');

corr_refract_length = -0.59;
dioptersPerMm = -1.24;

% Calculate the width of the distribution of axial lengths given the
% spherical refraction. This draws upon the conditional expectation of a
% X given Y in a bivariate normal distribution

sigma_length_givenRefraction = sqrt((1-corr_refract_length^2)*lengthGauss.c1);

% Create a set of 9 pupil radii which will be used for all the the searches
pupilRadii = 2+(randn(9,1)./5);

% Create variations in extrinsic translation vector and in axial length
resultIdx = 1;
for axialLengthDelta = -1:1:1
    for cameraX = -5:5:5
        for cameraY = -5:5:5
            for cameraZ = 125:25:175
                % We will create the ellipses using an axial length for the
                % eye that can is off by + and - 2SD of the conditional
                % distribution of axial lengths given knowledge of the
                % spherical refraction of the eye
                axialLength = 23.5924 + axialLengthDelta * 2 * sigma_length_givenRefraction;
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
%                estimatedSceneGeometry = estimateSceneGeometry(pupilData,'','useParallel',true,'ellipseArrayList',1:1:ellipseIdx-1);
                estimatedSceneGeometry='sd';
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


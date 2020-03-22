function [myObj, sceneGeometry] = createSceneObjective(videoStemName, frameSet, gazeTargets, setupArgs, keyVals)


%% Create initial sceneGeometry structure
sceneGeometry = createSceneGeometry(setupArgs{:});


%% Load the materials
load([videoStemName '_correctedPerimeter.mat'],'perimeter');
load([videoStemName '_glint.mat'],'glintData');
load([videoStemName '_pupil.mat'],'pupilData');


%% Assemble the args and key vals
% Extract the frames we want
perimeter.data = perimeter.data(frameSet);
glintData.X = glintData.X(frameSet); glintData.Y = glintData.Y(frameSet);
ellipseRMSE = pupilData.initial.ellipses.RMSE(frameSet);

% Assemble these components into the args variable
args = {perimeter, glintData, ellipseRMSE, gazeTargets};

%% Create the objective function
myObj = @(x) calcGlintGazeError( updateSceneGeometry( sceneGeometry, x ), args{:}, keyVals{:} );

%% Create a fixation function
% This anonymous function returns the [azi ele] of the eye for this
% sceneGeometry when the eye is fixated up

end





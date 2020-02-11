function [ objError, modelEyePose, modelGlint, modelPoseGaze, modelVecGaze, poseRegParams, vectorRegParams] = calcGlintGazeError( sceneGeometry, perimeter, gazeTargets, ellipseRMSE, glintData, modelEyePose, errorReg )
% The error in prediction of gaze and glint location for a sceneGeometry
%
% Syntax:
%   [ objError, gazeError, glintError, modelEyePose, modelGlint, modelGaze, regParams] = calcGlintGazeError( requestedError, sceneGeometry, perimeter, gazeTargets, ellipseRMSE, glintData, modelEyePose )
%
% Description:
%   The sceneGeometry defines a physical system that models the gaze
%   location of an eye based upon the perimeter of the pupil and the
%   location of a glint in the image. This routine returns the error in the
%   model prediction for these two components.
%
% Inputs:
%   requestedError        - String or char vector. Defines the type of
%                           error returned in the output variable objError.
%                           Options are {'gaze','glint','total'}
%   sceneGeometry
%   perimeter
%   gazeTargets
%   ellipseRMSE
%   glintData
%   modelEyePose          - Empty or 2xn vector, where n is the number of
%                           frames to be modeled. If left empty, the
%                           eyePose will be derived from the perimeter
%                           data.
%
% Optional key/value pairs:
%   none
%  'bar'                  - Scalar. Bar bar bar bar bar bar bar bar bar bar
%                           bar bar bar bar bar bar bar bar bar bar bar bar
%                           bar bar bar bar bar bar
%
% Outputs:
%   none
%   baz                   - Cell. Baz baz baz baz baz baz baz baz baz baz
%                           baz baz baz baz baz baz baz baz baz baz baz baz
%                           baz baz baz
%
% Examples:
%{
    load('/Users/aguirre/Dropbox (Aguirre-Brainard Lab)/TOME_processing/session2_spatialStimuli/TOME_3015/032417/EyeTracking/GazeCal03_correctedPerimeter.mat')
    load('/Users/aguirre/Dropbox (Aguirre-Brainard Lab)/TOME_processing/session2_spatialStimuli/TOME_3015/032417/EyeTracking/GazeCal03_glint.mat')
    load('/Users/aguirre/Dropbox (Aguirre-Brainard Lab)/TOME_processing/session2_spatialStimuli/TOME_3015/032417/EyeTracking/GazeCal03_pupil.mat')
    load('/Users/aguirre/Dropbox (Aguirre-Brainard Lab)/TOME_processing/session2_spatialStimuli/TOME_3015/032417/EyeTracking/GazeCal03_sceneGeometry.mat')

    frames = sceneGeometry.meta.estimateSceneParams.ellipseArrayList;
    perimeter.data = perimeter.data(frames);
	gazeTargets = sceneGeometry.meta.estimateSceneParams.fixationTargetArray;
    ellipseRMSE = pupilData.initial.ellipses.RMSE(frames);
    glintData.X = glintData.X(frames); glintData.Y = glintData.Y(frames);
    [ gazeError, glintError, modelEyePose, modelGlint, modelGaze] = ...
        calcGlintGazeError( sceneGeometry, perimeter, gazeTargets, ellipseRMSE, glintData )
    figure
    subplot(1,2,1);
    plot(gazeTargets(1,:),gazeTargets(2,:),'ok'); hold on;
    plot(modelGaze(1,:),modelGaze(2,:),'xr');
    xlabel('horizontal gaze [deg]'); ylabel('vertical gaze [deg]');
    title('Modeled gaze position');
    axis equal
    subplot(1,2,2);
    plot(glintData.X,glintData.Y,'ok'); hold on;
    plot(modelGlint.X,modelGlint.Y,'xr');
    xlabel('glint position [pixels]'); ylabel('glint position [pixels]');
    title('Modeled glint position');
    axis equal
%}


if nargin == 6
    errorReg = [2 1];
end

% How many frames do we have
nFrames = length(glintData.X);

% The weight for the errors is given by the inverse initial ellipse fit
% RMSE
weights = 1./ellipseRMSE;

% Handle if we are calculating eyePose
if isempty(modelEyePose)
    calcEyePose = true;
    modelEyePose = nan(nFrames,4);
else
    calcEyePose = false;
end

% Allocate the glint loop variables
modelGlintX = nan(nFrames,1);
modelGlintY = nan(nFrames,1);
perimFitError = nan(nFrames,1);
pupilCenter = nan(nFrames,2);
modelPupilCenter = nan(nFrames,2);

% These are some magic numbers used in retrieving the glint
nStopPerimPoints = 5;
glintIdx = nStopPerimPoints*2+1;

% Loop over the frames and obtain the modeled eyePose and glint
parfor ii = 1:nFrames
    
    % Get the perimeter
    Xp = perimeter.data{ii}.Xp;
    Yp = perimeter.data{ii}.Yp;
    
    % Get the glint
    glintCoord = [glintData.X(ii) glintData.Y(ii)];
    
    % Get the eyePose
    if calcEyePose
        [modelEyePose(ii,:),perimFitError(ii)] = eyePoseEllipseFit(Xp, Yp, sceneGeometry, 'glintCoord', glintCoord);
    end
        
    % Get the glint coordinates
    [fittedEllipse, fittedGlint] = projectModelEye(modelEyePose(ii,:), sceneGeometry, ...
        'nStopPerimPoints',nStopPerimPoints);
    
    % Store the empirical pupil center
    pupilCenter(ii,:) = fittedEllipse(1:2);
    
    % Store the glint coordinate, and Inf if we didn't get a glint
    if isempty(fittedGlint)
        modelGlintX(ii) = Inf;
        modelGlintY(ii) = Inf;
    else
        modelGlintX(ii) = fittedGlint(1);
        modelGlintY(ii) = fittedGlint(2);
    end
end

% Store the parpool loop variables
modelGlint.X = modelGlintX;
modelGlint.Y = modelGlintY;


%% Image error
% These are errors in matching the appearance of the eye in the image. The
% error is in pixel units

% perimError -- fit to the pupil perimeters
perimError = nanNorm(perimFitError,weights);

% glintError -- fit of the model to the glint locations
glintDistances = sqrt(sum([modelGlint.X - glintData.X, modelGlint.Y - glintData.Y].^2,2));
glintError = nanNorm(glintDistances,weights);

imageError = nanNorm([perimError glintError]);

%% Gaze error
% These are errors in matching the position of the fixation targets on the
% screen. The error is in degrees of visual angle. We need to have gaze
% targets to compute this error

if isempty(gazeTargets)
    poseError = nan;
    vectorError = nan;
    gazeError = nan;
else
    
    % poseError -- eye rotation equal to the visual angle
    poseRegParams = absor(...
        modelEyePose(:,1:2)',...
        gazeTargets,...
        'weights',weights,...
        'doScale',false,...
        'doTrans',true);
    modelPoseGaze = poseRegParams.R * modelEyePose(:,1:2)' + poseRegParams.t;
    poseError = nanNorm(sqrt(sum(gazeTargets - modelPoseGaze).^2)',weights);
    
    % vectorError -- vector between the glint and pupil center used to model
    % eye position
    glintSign = [1;-1];
    centerDiff = (pupilCenter - [modelGlint.X modelGlint.Y])' .* ...
        glintSign;
    if any(isnan(sum(centerDiff))) || any(isinf(sum(centerDiff)))
        vectorError = inf;
    else
        vectorRegParams = absor(...
            centerDiff,...
            gazeTargets,...
            'weights',weights,...
            'doScale',true,...
            'doTrans',true);
        vectorRegParams.glintSign = glintSign;
        modelVecGaze = vectorRegParams.s * vectorRegParams.R * centerDiff + vectorRegParams.t;
        vectorError = nanNorm(sqrt(sum(gazeTargets - modelVecGaze).^2)',weights);
    end
    
    gazeError = nanNorm([poseError, vectorError]);
end


%% Return the error
objError = nanNorm([imageError, gazeError],errorReg);
objError(isinf(objError))=realmax;

end

%%%%% LOCAL FUNCTIONS

function val = nanNorm(vec, weights)

% Prepare the weight vector
if nargin==1
    weights = ones(size(vec))./length(vec);
end
weights(isnan(vec))=nan;
weights = weights ./ nansum(weights);

% Calculate the normed value
val = sqrt(nansum((vec.^2).*weights));

end
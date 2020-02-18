function [ objError, modelEyePose, modelPupilEllipse, modelGlintCoord, modelPoseGaze, modelVecGaze, poseRegParams, vectorRegParams, rawErrors] = calcGlintGazeError( sceneGeometry, perimeter, glintData, ellipseRMSE, gazeTargets, varargin )
% Error in prediction of image and gaze for a sceneGeometry
%
% Syntax:
%   [ objError, modelEyePose, modelPupilEllipse, modelGlintCoord, modelPoseGaze, modelVecGaze, poseRegParams, vectorRegParams, rawErrors] = calcGlintGazeError( sceneGeometry, perimeter, glintData, ellipseRMSE, gazeTargets, varargin )
%
% Description:
%   The sceneGeometry defines a physical system of an eye, a camera (with a
%   light source), and fixation target. This routine returns the error in
%   the model prediction of  elements of this scene. Two of these
%   prediction elements are derived from the eye tracking image, and are
%   the error with which the perimeter of the pupil is fit with a
%   scene-constrained ellipse, and the error in the specification of the
%   location of the glint(s). Two more prediction elements relate to the
%   gaze location of the eye. If provided with a set of gaze targets in
%   degrees of visual angle, errors are calculate for matching the rotation
%   of the eye in degrees with the position of fixation targets, and for
%   matching the vector distance between the pupil center and glint to
%   these fixation targets.
%
% Inputs:
%   sceneGeometry         - Structure. SEE createSceneGeometry.m
%   perimeter             - Structure. SEE findPupilPerimeter.m
%   glintData             - Structure. SEE findPupilGlint.m
%   ellipseRMSE           - Vector. The error in the initial ellipse fit.
%   gazeTargets           - Vector.
%
% Optional key/value pairs:
%  'modelEyePose'         - Empty or 2xn vector, where n is the number of
%                           frames to be modeled. If left empty, the
%                           eyePose will be derived from the perimeter
%                           data. Passing the eyePose data speeds the
%                           execution of this routine.
%  'eyePoseLB/UB'         - 1x4 vector. Upper / lower bounds on the eyePose
%                           [azimuth, elevation, torsion, pupil radius].
%                           The torsion value is unusued and is bounded to
%                           zero. Biological limits in eye rotation and
%                           pupil size would suggest boundaries of [±35,
%                           ±25, 0, 0.25-4]. Note, however, that these
%                           angles are relative to the center of
%                           projection, not the primary position of the
%                           eye. Therefore, in circumstances in which the
%                           camera is viewing the eye from an off-center
%                           angle, the bounds will need to be shifted
%                           accordingly.
%  'errorReg'             - 1x4 vector. This vector defines a 
%                           regularization that weights the four types of
%                           error that can contribute to the overall model
%                           error. The weights apply to errors in:
%                               [perim glint pose vector]
%
% Outputs:
%   objError              - Scalar. The overall model error.
%   modelEyePose          - f x 4 matrix of modeled eyePose positions 
%                           for the passed frames.
%   modelPupilEllipse     - f x 5 matrix of ellipses fit to the pupil.
%   modelGlintCoord       - Structure. The modeled glint locations.
%   modelPoseGaze         - f x 2 matrix of modeled gaze locations in 
%                           horizontal and vertical degrees of visual angle
%                           derived from the eye rotation.
%   modelVecGaze          - f x 2 matrix of modeled gaze locations in 
%                           horizontal and vertical degrees of visual angle
%                           derived from the pupil center -> glint vec.
%   poseRegParams         - Structure. The parameters that relate eyePose
%                           to screen position.
%   vectorRegParams       - Structure. The parameters that relate the 
%                           pupil center -> glint vec to screen position.
%   rawErrors             - 1x4 matrix. The four component errors.
%
% Examples:
%{
%}



%% Parse input
p = inputParser;

% Required
p.addRequired('sceneGeometry',@isstruct);
p.addRequired('perimeter',@isstruct);
p.addRequired('glintData',@isstruct);
p.addRequired('ellipseRMSE',@isnumeric);
p.addRequired('gazeTargets',@isnumeric);

% Optional
p.addParameter('modelEyePose',[],@(x)(isempty(x) | isnumeric(x)));
p.addParameter('eyePoseLB',[-89,-89,0,0.1],@isnumeric);
p.addParameter('eyePoseUB',[89,89,0,4],@isnumeric);
p.addParameter('errorReg',[1 1 3 3],@isscalar);

% Parse and check the parameters
p.parse(sceneGeometry, perimeter, glintData, ellipseRMSE, gazeTargets, varargin{:});


%% Setup variables
% How many frames do we have
nFrames = length(glintData.X);

% The weight for the errors is given by the inverse initial ellipse fit
% RMSE
weights = 1./ellipseRMSE;

% Handle if we are calculating eyePose
if isempty(p.Results.modelEyePose)
    calcEyePose = true;
    modelEyePose = nan(nFrames,4);
else
    calcEyePose = false;
    modelEyePose = p.Results.modelEyePose;
end

% Allocate the loop variables
modelGlintX = nan(nFrames,1);
modelGlintY = nan(nFrames,1);
perimFitError = nan(nFrames,1);
pupilCenter = nan(nFrames,2);
modelPupilEllipse = nan(nFrames,5);

% These are some magic numbers used in retrieving the glint

% Loop over the frames and obtain the modeled eyePose and glint
parfor ii = 1:nFrames
    
    % Get the perimeter
    Xp = perimeter.data{ii}.Xp;
    Yp = perimeter.data{ii}.Yp;
    
    % Get the glint
    glintCoord = [glintData.X(ii) glintData.Y(ii)];
    
    % Get the eyePose
    if calcEyePose
        modelEyePose(ii,:) = eyePoseEllipseFit(Xp, Yp, sceneGeometry, 'glintCoord', glintCoord);
    end
        
    % Get the glint coordinates
    [modelPupilEllipse_loop, modelGlintCoord_loop] = projectModelEye(modelEyePose(ii,:), sceneGeometry);    
    modelPupilEllipse(ii,:) = modelPupilEllipse_loop;
    
    % Get the error in fitting the perimeter with the ellipse
    if any(isnan(modelPupilEllipse_loop))
        % Set fVal to something arbitrarily large
        perimFitError(ii) = 1e6;
    else
        % This is the RMSE of the distance values of the boundary
        % points to the ellipse fit.
        explicitEllipse = ellipse_transparent2ex(modelPupilEllipse_loop);
        if isempty(explicitEllipse)
            perimFitError(ii) = 1e6;
        else
            if any(isnan(explicitEllipse))
                perimFitError(ii) = 1e6;
            else
                perimFitError(ii) = sqrt(nanmean(ellipsefit_distance(Xp,Yp,explicitEllipse).^2));
            end
        end
    end
        
    % Store the empirical pupil center
    pupilCenter(ii,:) = modelPupilEllipse_loop(1:2);
    
    % Store the glint coordinate, and Inf if we didn't get a glint
    if isempty(modelGlintCoord_loop)
        modelGlintX(ii) = Inf;
        modelGlintY(ii) = Inf;
    else
        modelGlintX(ii) = modelGlintCoord_loop(1);
        modelGlintY(ii) = modelGlintCoord_loop(2);
    end
end

% Store the parpool loop variables
modelGlintCoord.X = modelGlintX;
modelGlintCoord.Y = modelGlintY;


%% Image error
% These are errors in matching the appearance of the eye in the image. The
% error is in pixel units

% perimError -- fit to the pupil perimeters
perimError = nanNorm(perimFitError,weights);

% glintError -- fit of the model to the glint locations
glintDistances = sqrt(sum([modelGlintCoord.X - glintData.X, modelGlintCoord.Y - glintData.Y].^2,2));
glintError = nanNorm(glintDistances,weights);


%% Gaze error
% These are errors in matching the position of the fixation targets on the
% screen. The error is in degrees of visual angle. We need to have gaze
% targets to compute this error

if isempty(gazeTargets)
    poseError = nan;
    vectorError = nan;
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
    centerDiff = (pupilCenter - [modelGlintCoord.X modelGlintCoord.Y])' .* ...
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
    
end


%% Return the error
rawErrors = [perimError glintError poseError, vectorError];
objError = nanNorm(rawErrors,p.Results.errorReg);
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
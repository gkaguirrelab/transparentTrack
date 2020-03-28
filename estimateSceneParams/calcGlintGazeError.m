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
%  'missedGlintPenalty'   - Scalar. If a glint is not found for one frame,
%                           then this value is assigned as the error for
%                           that frame.
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
p.addParameter('eyePoseLB',[-89,-89,0,0.1],@isnumeric);
p.addParameter('eyePoseUB',[89,89,0,4],@isnumeric);
p.addParameter('errorReg',[1 2 4 2],@isnumeric);
p.addParameter('missedGlintPenalty',1e3,@isnumeric);
p.addParameter('poseRegParams',[],@isstruct);
p.addParameter('vectorRegParams',[],@isstruct);


% Parse and check the parameters
p.parse(sceneGeometry, perimeter, glintData, ellipseRMSE, gazeTargets, varargin{:});


%% Setup variables
% How many frames do we have
nFrames = length(glintData.X);

% The weight for the errors is given by the inverse initial ellipse fit
% RMSE
weights = 1./ellipseRMSE;


% Allocate the loop and return variables
modelEyePose = nan(nFrames,4);
modelGlintX = nan(nFrames,1);
modelGlintY = nan(nFrames,1);
perimFitError = nan(nFrames,1);
pupilCenter = nan(nFrames,2);
modelPupilEllipse = nan(nFrames,5);
modelPoseGaze = nan(nFrames,2);
modelVecGaze = nan(nFrames,2);
poseRegParams = struct();
vectorRegParams = struct();


% Loop over the frames and obtain the modeled eyePose and glint
parfor ii = 1:nFrames
    
    % Get the perimeter
    Xp = perimeter.data{ii}.Xp;
    Yp = perimeter.data{ii}.Yp;
    
    % Get the glint
    glintCoord = [glintData.X(ii) glintData.Y(ii)];
    
    % Get the eyePose
    modelEyePose(ii,:) = eyePoseEllipseFit(Xp, Yp, sceneGeometry, 'glintCoord', glintCoord);
        
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
glintDistances(isinf(glintDistances)) = p.Results.missedGlintPenalty;
glintError = nanNorm(glintDistances,weights);


%% Gaze error
% These are errors in matching the position of the fixation targets on the
% screen. The error is in degrees of visual angle. We need to have gaze
% targets to compute this error
if isempty(gazeTargets)
    poseError = nan;
    vectorError = nan;
else
    %% poseError
    % Match eye rotation to visual angle of fixation targets
    
    % Identify frames for which we have a gaze target, and thus are
    % eligible for this calculation
    validFrames = ~isnan(sum(gazeTargets));
    
    % Create variable to hold the error per frame and the modeled gaze. By
    % default, a frame that does not have an associate gaze target will not
    % contribute to the calculation of the poseError, so this vector is
    % initialized as nans.
    poseErrorByFrame = nan(size(sum(gazeTargets)));
    modelPoseGaze = nan(size(modelEyePose(:,1:2)))';
    
    % Either calculate or use the supplied poseRegParams
    if isempty(p.Results.poseRegParams)
        poseRegParams = absor(...
            modelEyePose(validFrames,1:2)',...
            gazeTargets(:,validFrames),...
            'weights',weights(validFrames),...
            'doScale',false,...
            'doTrans',true);
    else
        poseRegParams = p.Results.poseRegParams;
    end
    
    % Calculate the eyePose gaze and error for the modeled frames
    modelPoseGaze(:,validFrames) = poseRegParams.R * modelEyePose(validFrames,1:2)' + poseRegParams.t;
    poseErrorByFrame(validFrames) = sqrt(sum( (gazeTargets(:,validFrames) - modelPoseGaze(:,validFrames)).^2 ))';
    
    % Store the pose error
    poseError = nanNorm(poseErrorByFrame',weights);
    
    
    %% vectorError
    % Find a transformation between gaze pose and the vectors that connect
    % the glint and pupil center
    
    % The relationship between image y glint position and eye elevation is
    % sign reversed. Thinking this through is left as an exercise to the
    % reader.
    glintSign = [1;-1];
    
    % Obtain the x,y differences between the pupil center and glint for
    % each frame
    centerDiff = (pupilCenter - [modelGlintCoord.X modelGlintCoord.Y])' .* ...
        glintSign;
    
    % Valid frames are those for which we have been supplied a gaze target,
    % and those for which the projection model is able to calculate a
    % predicted glint location. Those frames for which we lack a gazePose
    % are marked nan and excluded from the error calculation. 
    % Defin

    % Create variable to hold the error per frame and the modeled gaze. By
    % default, a frame that does not have an associate gaze target will not
    % contribute to the calculation of the poseError, so this vector is
    % initialized as nans. 
    vectorErrorByFrame = nan(size(sum(gazeTargets)));
    modelVecGaze = nan(size(modelEyePose(:,1:2)))';    
    
    % Those frames for which the model cannot produce a predicted glint
    % contribute to the objective function error. If the model cannot
    % produce a predicted glint, then the scene parameters are extreme /
    % invalid. Therefore, we assign the missedGlintPenalty to those frames
    % for which the model cannot calculate a predicted glint.
    noGlintFrames = ~isfinite(sum(centerDiff));
    vectorErrorByFrame(noGlintFrames) = p.Results.missedGlintPenalty;
    
    % The gaze position is calculated for the valid frames, which are those
    % frames for which we have both a gaze target and a predicted glint
    % produced by the projection model.
    validFrames = and(~isnan(sum(gazeTargets)),~noGlintFrames);

    % Either calculate or use the supplied vectorRegParams
    if isempty(p.Results.vectorRegParams)
        vectorRegParams = absor(...
            centerDiff(:,validFrames),...
            gazeTargets(:,validFrames),...
            'weights',weights(validFrames),...
            'doScale',true,...
            'doTrans',true);
        vectorRegParams.glintSign = glintSign;
    else
        vectorRegParams = p.Results.vectorRegParams;
    end
        
    % Calculate the vector-based gaze and error for the modeled frames
    modelVecGaze(:,validFrames) = vectorRegParams.s * vectorRegParams.R * centerDiff(:,validFrames) + vectorRegParams.t;
    vectorErrorByFrame(validFrames) = sqrt(sum( (gazeTargets(:,validFrames) - modelVecGaze(:,validFrames)).^2 ));

    % Store the vector error
    vectorError = nanNorm(vectorErrorByFrame',weights);
end


%% Return the error
rawErrors = [perimError glintError poseError, vectorError];
objError = nanNorm(rawErrors,p.Results.errorReg);
objError(isinf(objError))=realmax;


end


%%%%% LOCAL FUNCTIONS

function val = nanNorm(vec, weights)
% Returns a weighted Euclidean norm, ignoring nan values
%
% Syntax:
%  val = nanNorm(vec, weights)
%

% Prepare the weight vector
if nargin==1
    weights = ones(size(vec))./length(vec);
end
weights(isnan(vec))=nan;
weights = weights ./ nansum(weights);

% Calculate the normed value
val = sqrt(nansum((vec.^2).*weights));

end
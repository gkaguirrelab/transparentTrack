function updateError( obj, varargin )
% Calculate the error in prediction of image and gaze for the sceneGeometry
%
% Syntax:
%  obj.updateError()
%
% Description:
%   The sceneGeometry defines a physical system of an eye, a camera (with a
%   light source), and fixation target. This routine calculates the error
%   in the model prediction of elements of this scene. Two of these
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
%   None. These are obtained from the obj.
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
%  'poseRegParams', 'vectorRegParams' - Structures. These values, generated
%                           from a previous execution of this function, can
%                           be returned in the varargin to be used instead
%                           of being recomputed. This is done to constrain
%                           the model, as opposed to saving computation
%                           time.
%
% Outputs:
%   None, but these elements are updated in the obj:
%
%   fVal                  - Scalar. The overall model error.
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


%% Obtain variables from the object
sceneGeometry = obj.sceneGeometry;
perimeter = obj.perimeter;
glintDataX = obj.glintDataX;
glintDataY = obj.glintDataY;
ellipseRMSE = obj.ellipseRMSE;
gazeTargets = obj.gazeTargets;
relCamPos = obj.relCamPos(:,obj.frameSet);


%% Parse input
p = inputParser;

% Optional
p.addParameter('eyePoseLB',[-89,-89,0,0.1],@isnumeric);
p.addParameter('eyePoseUB',[89,89,0,4],@isnumeric);
p.addParameter('errorReg',[1 2 4 2],@isnumeric);
p.addParameter('missedGlintPenalty',1e3,@isnumeric);
p.addParameter('poseRegParams',[],@isstruct);
p.addParameter('vectorRegParams',[],@isstruct);

% Parse and check the parameters
p.parse(varargin{:});


%% Setup variables

% How many frames do we have
nFrames = length(glintDataX);

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


%% Loop over the frames and obtain the modeled eyePose and glint
parfor ii = 1:nFrames
    
    % Get the perimeter
    Xp = perimeter{ii}.Xp;
    Yp = perimeter{ii}.Yp;
    
    % Get the glint
    glintCoord = [glintDataX(ii) glintDataY(ii)];
    
    % Update the sceneGeometry camera position to account for known head
    % translation
    adjustedSceneGeometry = sceneGeometry;
    cameraPosition = sceneGeometry.cameraPosition.translation;
    cameraPosition = cameraPosition + relCamPos(:,ii);
    adjustedSceneGeometry.cameraPosition.translation = cameraPosition;
    
    % Get the eyePose
    modelEyePose(ii,:) = eyePoseEllipseFit(Xp, Yp, adjustedSceneGeometry, 'glintCoord', glintCoord);
    
    % Get the glint coordinates
    [modelPupilEllipse_loop, modelGlintCoord_loop] = projectModelEye(modelEyePose(ii,:), adjustedSceneGeometry);
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
glintDistances = sqrt(sum([modelGlintCoord.X - glintDataX, modelGlintCoord.Y - glintDataY].^2,2));
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

        % Add a meta field to the params with the formula
        poseRegParams.meta = 'f = R * [azi; ele] + t';

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
        
        % Add a meta field to the params with the formula
        vectorRegParams.meta = 'f = s * R * [Xp - Xg; -1*(Yp-Yg)] + t';
    else
        vectorRegParams = p.Results.vectorRegParams;
    end
    
    % Calculate the vector-based gaze and error for the modeled frames
    modelVecGaze(:,validFrames) = vectorRegParams.s * vectorRegParams.R * centerDiff(:,validFrames) + vectorRegParams.t;
    vectorErrorByFrame(validFrames) = sqrt(sum( (gazeTargets(:,validFrames) - modelVecGaze(:,validFrames)).^2 ));
    
    % Store the vector error
    vectorError = nanNorm(vectorErrorByFrame',weights);
end


%% Obtain the omnibus error
rawErrors = [perimError glintError poseError, vectorError];
fVal = nanNorm(rawErrors,p.Results.errorReg);
fVal(isinf(fVal))=realmax;


%% Store the results in the object

% Store the fVal
obj.fVal = fVal;

% Store all the other model components
obj.modelEyePose = modelEyePose;
obj.modelPupilEllipse = modelPupilEllipse;
obj.modelGlintCoord = modelGlintCoord;
obj.modelPoseGaze = modelPoseGaze;
obj.modelVecGaze = modelVecGaze;
obj.poseRegParams = poseRegParams;
obj.vectorRegParams = vectorRegParams;
obj.rawErrors = rawErrors;

% Save the fixation
obj.fixationEyePose = poseRegParams.t;
obj.screenTorsion = poseRegParams.theta;
obj.screenRotMat = poseRegParams.R;

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
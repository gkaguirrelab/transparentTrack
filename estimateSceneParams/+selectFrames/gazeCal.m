function [frameSet, gazeTargets, eyePose, rho, theta] = gazeCal(videoStemName)
% Identify the fixation frame that can be used to sync sceneGeometry
%
% Syntax:
%  [frameSet, gazeTargets] = selectFrames.gazeCal(videoStemName)
%
% Description:
%   Positioning an eye model in a scene requires the selection of
%   informative frames of the acquisition to guide the alignment.
%
%   This function identifies the frame from a gazeCal acquisition that
%   corresponds to fixation at the [0;0] position.
%
% Inputs:
%	videoStemName         - Char vector. Full path to video file from which
%                           the scene observations have been derived. The
%                           stem name should omit the "_gray.avi" suffix
%                           that is usually present in the names of these
%                           video files.
%
% Optional key-value pairs:
%   none
%
% Outputs:
%   frameSet              - Scalar that specifies a frame index
%                           (indexed from 1).
%   gazeTargets           - A 2x1 matrix that provides the positions, in
%                           degrees of visual angle, of the likely fixation
%                           position [0;0] of the eye for this frame.
%



%% input parser
p = inputParser; p.KeepUnmatched = true;

% Required
p.addRequired('videoStemName',@ischar);

% parse
p.parse(videoStemName)


% Load the sceneGeometry that already exists for this gazeCal acquisition
load([videoStemName '_sceneGeometry.mat'],'sceneGeometry');

% Get the frameSet and gazeTargets
frameSet = sceneGeometry.meta.estimateSceneParams.obj.frameSet;
gazeTargets = sceneGeometry.meta.estimateSceneParams.obj.gazeTargets;

% Which of the list of frames is the [0;0] fixation frame?
idx = logical((gazeTargets(1,:)==0).*(gazeTargets(2,:)==0));

% Define these in case we return early
rho = [];
theta = [];

% Some gazeCal runs lacked a formal measurement of the [0;0] target. Check
% if we have one.
if sum(idx)==1

    % We do. Return the fixation frame
    frameSet = frameSet(idx);
    gazeTargets = gazeTargets(:,idx);
    
    % Get the modeled eyePose for this frame
    eyePose = sceneGeometry.meta.estimateSceneParams.obj.modelEyePose(idx,:);

    % Get the pupil ellipse
    pupilEllipseFixationIn = sceneGeometry.meta.estimateSceneParams.obj.modelPupilEllipse(idx,:);

    % Find the shape of the pupil for this frame, expressed as theta and rho
    % values (SEE: csaEllipseError)
    rho = 1-sqrt(1-pupilEllipseFixationIn(4)^2);
    theta = pupilEllipseFixationIn(5)*2;
    
else
    % We don't. Let's see what else we can do...
    
    % Calculate the eyePose that corresponds to the predicted fixation of
    % the [0;0] screen position.
    R = sceneGeometry.screenPosition.poseRegParams.R;
    t = sceneGeometry.screenPosition.poseRegParams.t;
    g = R\(-t);    
    eyePose = [g(1) g(2) 0 2];

    % See if we have eyePoses calculated for the gazeCal acquisition. If
    % so, we can find the frame with the eyePose closest to g
    pupilFileName = [videoStemName '_pupil.mat'];
    if isfile(pupilFileName)
        load(pupilFileName,'pupilData');
        if isfield(pupilData,'sceneConstrained')
            % Find the frame with the closest eyePose            
            [~,frameSet] = min(vecnorm(g - pupilData.sceneConstrained.eyePoses.values(:,1:2)'));
            gazeTargets = [0;0];

            % Obtain the rho and theta values of the pupil ellipse for this gaze
            % position
            pupilEllipse = projectModelEye(pupilData.sceneConstrained.eyePoses.values(frameSet,:),sceneGeometry);
            rho = pupilEllipse(4);
            theta = pupilEllipse(5);
            return
        end
    end

    % If we are still in this function, we will have to select a frame
    % based upon the shape of the pupil ellipse.
    
    % Obtain the rho and theta values of the pupil ellipse for this gaze
    % position
    pupilEllipse = projectModelEye(eyePose,sceneGeometry);    
    rho = pupilEllipse(4);
    theta = pupilEllipse(5);
    
    % Find the frame with this pupil shape 
    [frameSet, gazeTargets] = selectFrames.shape(videoStemName, rho, theta);

    % Issue a warning that we could not find a true fixation frame
    warning('transparentTrack:selectFrames_gazeCal','No fixation frame for this gazeCal; consider creating a sceneConstrained pupil file');
        
end

end
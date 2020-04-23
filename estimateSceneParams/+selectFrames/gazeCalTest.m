function [frameSet, gazeTargets] = gazeCalTest(videoStemName, varargin)
% Identify a fixation frame that can be used to sync sceneGeometry
%
% Syntax:
%  [frameSet, gazeTargets] = selectFrames.gazePre(videoStemName)
%
% Description:
%   Positioning an eye model in a scene requires the selection of
%   informative frames of the acquisition to guide the alignment.
%
%   This function is used for the particular circumstance in which we are
%   testing the behavior of the sceneGeometry routines by syncing one
%   gazeCal acquisition to another. This function identifies the frame in a
%   the acquisition that corresponds to fixation of the [0 0] fixation
%   target.
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
p.parse(videoStemName, varargin{:})


% Load the sceneGeometry that already exists for this gazeCal acquisition
load([videoStemName '_sceneGeometry.mat'],'sceneGeometry');

% Which of the list of frames is the [0;0] fixation frame?
idx = find((sceneGeometry.meta.estimateSceneParams.obj.gazeTargets(1,:)==0).*(sceneGeometry.meta.estimateSceneParams.obj.gazeTargets(2,:)==0));

% Some gazeCal runs lacked a formal measurement of the [0;0] target. Check
% if we have one.
if ~isempty(idx)

    % We do. Identify the absolute index of the frame
    frameSet = sceneGeometry.meta.estimateSceneParams.obj.frameSet(idx);

else
    % We don't. In this case, grab a frame based upon pupil
    % shape predicted by the sceneGeometry for the [0;0] position
    
    % Calculate the eyePose that corresponds to the predicted fixation of
    % the [0;0] screen position.
    R = sceneGeometry.screenPosition.poseRegParams.R;
    t = sceneGeometry.screenPosition.poseRegParams.t;
    g = inv(R)*(-t);
    
    % Obtain the rho and theta values of the pupil ellipse for this gaze
    % position
    eyePose = [g(1) g(2) 0 2];
    pupilEllipse = projectModelEye(eyePose,sceneGeometry);    
    rho = pupilEllipse(4);
    theta = pupilEllipse(5);
    
    % Find the frame with this pupil shape 
    frameSet = shape(videoStemName, rho, theta);
    
end

% The gaze target is presumed to be the center of the screen [0 0]
gazeTargets = [0; 0];

end
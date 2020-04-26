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
%   gazeCal acquisition to another. This function identifies the frames in
%   the acquisition that are the gaze fixation frames, and returns these
%   sorted so that the frame with fixation at the [0;0] position is first.
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


% This routine identifies a fixation frame
[frameSetFix, gazeTargetsFix] = select.gazeCal(videoStemName);

% Load the sceneGeometry that already exists for this gazeCal acquisition
load([videoStemName '_sceneGeometry.mat'],'sceneGeometry');

% Get the frameSet and gazeTargets
frameSet = sceneGeometry.meta.estimateSceneParams.obj.frameSet;
gazeTargets = sceneGeometry.meta.estimateSceneParams.obj.gazeTargets;

% Which of the list of frames is the [0;0] fixation frame?
idx = logical((gazeTargets(1,:)==0).*(gazeTargets(2,:)==0));

% Sort the return variables so that the fixation frame is first
frameSet = [frameSetFix, frameSet(~idx)];
gazeTargets = [gazeTargetsFix, gazeTargets(:,~idx)];


end
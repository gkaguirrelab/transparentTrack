function calcCalibrationMatrix (gazeDataFileName,calMatrixFileName,varargin)

% calcCalibrationMatrix (gazeDataFileName,calMatrixFileName)
% 
% this function calculates the gaze calibration matrix starting from a set
% of known targets observed at a known distance, and apparent gaze
% direction for each fixation. The calibration matrix is the one that
% minimizes the distance between all the target locations and the
% projection of the apparent gaze direction for each fixation.
% 
% The calibration matrix is a 4x4 matrix in homogeneous coordinates and can
% be applied to the raw pupil and glint data as follows:
% 
% [aXYZW] = calMatrix * [(pX-gX)/Rpc; (pY-gY)/Rpc; (1 - sqrt(((pX-gX)/Rpc)^2 + ((pY-gY)/Rpc)^2)); 1];
% [calGazeX;calGazeY;viewingDistance] =    (aXYZW(1:3)/aXYZW(4))';
% 
% Where:
%   aXYZW = calibrated Gaze Data in homogeneous screen coordinates
%   [pX pY] = center of pupil in pixels
%   [gX gY] = center of glint in pixels
%   Rpc = relative perspective correction (between target and apparent gaze
%       vector)
%   
% Note that the first line applies the calMatrix to the 3-D projection of
% the apparent gaze vector (in pixels) in homogeneous coordinates, while
% the second line converts the calibrated data from homogeneous
% coordinates to 3-D screen coordinates, where the 3rd dimension is the
% distance between the observer and the screen.
%
% 
% MORE ABOUT CALCULATING THE CALIBRATION MATRIX
% 
% The calibration matrix converts the apparent gaze vector (center of pupil
% - center of glint) from pixel coordinates (origin top left corner of the
% frame) to "screen coordinates".
% 
% Screen coordinates are defined as millimeters from the center of the
% screen (X growing left to right, Y growing top to bottom) when the
% subject is looking at the screen from a defined viewing distance (also in
% mm).
% 
% The transformation from 3-D screen coordinates to 2-D pixel coordinates
% can be written in homogeneous coordinates in the most general form as
% follows:
% 
% [px; py; 1] = M * [Sx; Sy; viewDist; 1]
% 
%  where M is a 3x4 transformation matrix operating: rotation, scaling,
%  stretching and 2-d projection on the screen coordinates.
% 
% In this case we break the conversion problem in the following steps:
% 
% 1. estimate the relative perspective correction factor (Rpc) between the
% appearance of the targets and of the apparent gaze vectors. Note that
% this factor has nothing to do with the geometry of the eye! It solely
% depends on the difference between the distance of targets and apparent
% gaze from the camera plane.
% 
% 2. define the 3-D projection in pixel coordinates of the apparent gaze
% vector, using an estimate of the distance from center of corneal
% curvature to the pupil in pixel in the definition of the 3rd dimension.
% 
% 3. express both target and apparent gaze in homogeneous coordinates.
% 
% 4. Correct the apparent gaze vector in homogeneous coordinates applying a
% perspective division using Rpc as divisor (this step is also called
% homogeneous divide). More info on perspective correction in the header of
% calcRpc.
% 
% 5. use fminsearch to find the matrix that minimizes the total distance
% between the target location and the apparent gaze location for each
% fixation.


%% input parser

p = inputParser; p.KeepUnmatched = true;

% Required
p.addRequired('gazeDataFileName',@ischar);
p.addRequired('calMatrixFileName',@ischar);

% Optional analysis parameters


% Optional display and I/O parameters
p.addParameter('verbosity','none', @ischar);

% Environment parameters
p.addParameter('tbSnapshot',[],@(x)(isempty(x) | isstruct(x)));
p.addParameter('timestamp',char(datetime('now')),@ischar);
p.addParameter('username',char(java.lang.System.getProperty('user.name')),@ischar);
p.addParameter('hostname',char(java.net.InetAddress.getLocalHost.getHostName),@ischar);

% parse
p.parse(gazeDataFileName, calMatrixFileName, varargin{:})


%% load gaze data
gazeData = load(gazeDataFileName);

targets = [gazeData.gazeCalData.targets.X, gazeData.gazeCalData.targets.Y];
pupil = [gazeData.gazeCalData.pupil.X, gazeData.gazeCalData.pupil.Y];
glint = [gazeData.gazeCalData.glint.X, gazeData.gazeCalData.glint.Y];
viewingDistance = gazeData.gazeCalData.viewingDistance;

clear gazeData

%% Estimate the relative perspective correction factor



% theoretically, the Z term of the apparent gaze should include the
% distance between the pupil center and the center of the corneal
% curvature in pixels (is ranges between 50-200). Given the subsequent steps of dividing by the Rpc and of
% fminsearch, we can set distancePupilCornea/Rpc as 1 without loss of
% generality.



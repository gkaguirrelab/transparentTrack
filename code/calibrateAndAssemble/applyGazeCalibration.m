function [calibratedGaze] = applyGazeCalibration(pupilFileName,glintFileName,gazeCalFactorsFileName,varargin)
% Applies the gaze calibration parameters to the raw pupil
% and glint data
% 
% Description:
%    This function applies the gaze calibration parameters to the raw pupil
%    and glint data, using the "gaze vector" approach. The gaze vector is
%    defined as (P - G) , where P and G are the X,Y position on the scene of 
%    the center of the pupil and the glint, respectively.
%   
%    All the necessary calibration parameters are stored in the gazeCalParams
%    file, and are derived from one or more trials in which the gaze vector is
%    calculated for targets presented on screen at known locations.
%    
%    The calibration strategy is as follows:
% 
%    [aXYZW] = calMatrix * [(pX-gX)/Rpc; (pY-gY)/Rpc; (1 - sqrt(((pX-gX)/Rpc)^2 + ((pY-gY)/Rpc)^2)); 1];
%    [calGazeX;calGazeY;viewingDistance] =    (aXYZW(1:3)/aXYZW(4))';
% 
%    Where:
%      aXYZW = calibrated Gaze Data in homogeneous screen coordinates
%      calMatrix = calibration matrix
%      [pX pY] = center of pupil in pixels
%      [gX gY] = center of glint in pixels
%      Rpc = relative perspective correction (between target and apparent gaze
%            vector)
%   
% 	 Note that the first line applies the calMatrix to the 3-D projection of
% 	 the apparent gaze vector (in pixels) in homogeneous coordinates, while
% 	 the second line converts the calibrated data from homogeneous
% 	 coordinates to 3-D screen coordinates, where the 3rd dimension is the
%    distance between the observer and the screen.
% 
% 	 More info about the calibration strategy can be found in the
%    calcCalibrationMatrix and calcRpc functions header.
% 
% Inputs:
%   pupilFileName           - name of the file with the pupil data to be calibrated, 
%                             as it results from the pupil pipeline.
%   glintFileName           - name of the mat with glint data
%   gazeCalFactorsFileName  - name of the mat file with the gaze calibration
%                             params.
% 
% Optional key/value pairs (display and I/O)
%  'verbosity' 			`   - level of verbosity. [none, full]
%
% Optional key/value pairs (environment)
%  'tbSnapshot'             - This should contain the output of the
%                             tbDeploymentSnapshot performed upon the result
%                             of the tbUse command. This documents the state
%                             of the system at the time of analysis.
%  'timestamp'              - AUTOMATIC; The current time and date
%  'username'               - AUTOMATIC; The user
%  'hostname'               - AUTOMATIC; The host
%
% Optional key/value pairs (analysis)
%  'calibratedGazeFileName' - name of the output file containing the
%       					  calibrated data, if the user wishes to save it on file.
%   'calibratedUnits'       - units in which the calibrated data is expressed
%                             (default [mmOnScreen])
%   'whichFitToCalibrate'   - which of the pupil fit in the pupil file to calibrate
%   	                     (default 'radiusSmoothed', also available 'sceneConstrained' and 'initial')
%   'analysisPass'          - set the pass number in case calibration data undergoes
%                             some kind of iterative correction process.
%
%   Outputs:
%   calibratedGaze          - structure containing the calibrated pupil width, height
%                             and area. The calibrated units are dependent on the size calibration
%                    	      method used.
% 

%% Parse vargin for options passed here
p = inputParser; p.KeepUnmatched = true;

% Required
p.addRequired('pupilFileName',@ischar);
p.addRequired('glintFileName',@ischar);
p.addRequired('gazeCalFactorsFileName',@ischar);

% Optional analysis parameters
p.addParameter('calibratedGazeFileName','',@ischar);
p.addParameter('calibratedUnits','mm', @ischar);
p.addParameter('whichFitToCalibrate','radiusSmoothed', @ischar);
p.addParameter('analysisPass',1, @isnumeric);

% Optional display and I/O parameters
p.addParameter('verbosity','none', @ischar);

% Environment parameters
p.addParameter('tbSnapshot',[],@(x)(isempty(x) | isstruct(x)));
p.addParameter('timestamp',char(datetime('now')),@ischar);
p.addParameter('username',char(java.lang.System.getProperty('user.name')),@ischar);
p.addParameter('hostname',char(java.net.InetAddress.getLocalHost.getHostName),@ischar);

% parse
p.parse(pupilFileName, glintFileName, gazeCalFactorsFileName, varargin{:})

%% load pupil data
tmpData = load(pupilFileName);
% pull transparent raw pupil data
rawPupilTransparent = tmpData.pupilData.(p.Results.whichFitToCalibrate);
% extract X and Y coordinates of pupil center
pupil.X = rawPupilTransparent.ellipse(:,1);
pupil.Y = rawPupilTransparent.ellipse(:,2);

clear tmpData
clear rawPupilTransparent
%% load glint data
tmpData = load(glintFileName);
% extract X and Y coordinates of glint center
glint.X = tmpData.glintData.X;
glint.Y = tmpData.glintData.Y;

clear tmpData
%% load gaze calibration params
tmpData = load(gazeCalFactorsFileName);

calMatrix = tmpData.gazeCalFactors.calMatrix;
perspectiveCorrection = tmpData.gazeCalFactors.perspectiveCorrection;

clear tmpData
%% apply calibration to data

% intialize gaze variables
calibratedGaze.X = nan(size(pupil.X));
calibratedGaze.Y = nan(size(pupil.X));
screenCoords = nan(length(pupil.X),3);
polarCoords = nan(length(pupil.X),2);

for ii = 1:length(pupil.X)
    pX = pupil.X(ii);
    pY = pupil.Y(ii);
    gX = glint.X(ii);
    gY = glint.Y(ii);
    % calculate apparent gaze vector in homogeneous coordinates
    aXYZW = calMatrix * [...
        (pX-gX)/perspectiveCorrection; ...
        (pY-gY)/perspectiveCorrection; ...
        (1 - sqrt(((pX-gX)/perspectiveCorrection)^2 + ((pY-gY)/perspectiveCorrection)^2)); ...
        1];
    % homogeneous divide (convert in 3-d screen coordinates)
    screenCoords(ii,:) = (aXYZW(1:3)/aXYZW(4))';
    
    % convert in polar coordinates (with same viewing distance)
    [polarCoords(ii,1), polarCoords(ii,2)] = convertScreenCoordinatesToPolar(screenCoords(ii,1),screenCoords(ii,2),screenCoords(ii,3));
end

%  pull out the coordinates
calibratedGaze.X = screenCoords(:,1);
calibratedGaze.Y = screenCoords(:,2);
calibratedGaze.ecc = polarCoords(:,1);
calibratedGaze.pol = polarCoords(:,2);
calibratedGaze.viewingDist = screenCoords(:,3);

%% save out calibrated gaze and metadata
calibratedGaze.meta = p.Results;

if ~isempty(p.Results.calibratedGazeFileName)
    save(p.Results.calibratedGazeFileName,'calibratedGaze');
end


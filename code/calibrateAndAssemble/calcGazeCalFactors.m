function calcGazeCalFactors (gazeDataFileName,gazeCalFactorsFileName,varargin)
% calculates the gaze calibration matrix starting from a set of known targets
% 
% Description:
%   this function calculates the gaze calibration matrix starting from a
%   set of known targets observed at a known distance, and apparent gaze
%   direction for each fixation. The calibration matrix is the one that
%   minimizes the distance between all the target locations and the
%   projection of the apparent gaze direction for each fixation.
% 
%   The calibration matrix is a 4x4 matrix in homogeneous coordinates and
%   can be applied to the raw pupil and glint data as follows:
% 
% [aXYZW] = calMatrix * [(pX-gX)/perspectiveCorrection; (pY-gY)/perspectiveCorrection; (1 - sqrt(((pX-gX)/perspectiveCorrection)^2 + ((pY-gY)/perspectiveCorrection)^2)); 1];
% [calGazeX;calGazeY;viewingDistance] =    (aXYZW(1:3)/aXYZW(4))';
% 
%   Where:
%       aXYZW = calibrated Gaze Data in homogeneous screen coordinates
%       [pX pY] = center of pupil in pixels
%       [gX gY] = center of glint in pixels
%       perspectiveCorrection = perspective correction (between target and apparent gaze
%           vector)
%   
%   Note that the first line applies the calMatrix to the 3-D projection of
%   the apparent gaze vector (in pixels) in homogeneous coordinates, while
%   the second line converts the calibrated data from homogeneous
%   coordinates to 3-D screen coordinates, where the 3rd dimension is the
%   distance between the observer and the screen.
%
% Input (required)
%   gazeDataFileName       - name of the mat file containing the gaze calibration
%                            data
%   gazeCalFactorsFileName - name of the mat file to save the calibration
%                           factors
% 
% Optional key/value pairs (analysis)
%   'fminsearchCalls'      - number of iteration for fminseach, changing the
%                            tolerance
%   'showFigures'          - toggle on to show the figure with the results
%                            of the calibration test
%
% Optional key/value pairs (display and I/O)
%  'verbosity'             - level of verbosity. [none, full]
%
% Optional key/value pairs (environment)
%   'tbSnapshot'           - This should contain the output of the
%                            tbDeploymentSnapshot performed upon the result
%                            of the tbUse command. This documents the state
%                            of the system at the time of analysis.
%   'timestamp'            - AUTOMATIC; The current time and date
%   'username'             - AUTOMATIC; The user
%   'hostname'             - AUTOMATIC; The host
% 
% Output (saved to file)
%   gazeCalFactors         - struct containing the calibration matrix and
%                            the perspective correction value to calibrate
%                            raw data as shown above.
% 
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
% 1. estimate the relative perspective correction factor
% (perspectiveCorrection) between the appearance of the targets and of the
% apparent gaze vectors. Note that this factor has nothing to do with the
% geometry of the eye! It solely depends on the difference between the
% distance of targets and apparent gaze from the camera plane.
% 
% 2. define the 3-D projection in pixel coordinates of the apparent gaze
% vector, using an estimate of the distance from center of corneal
% curvature to the pupil in pixel in the definition of the 3rd dimension.
% 
% 3. express both target and apparent gaze in homogeneous coordinates.
% 
% 4. Correct the apparent gaze vector in homogeneous coordinates applying a
% perspective division using perspectiveCorrection as divisor (this step is also called
% homogeneous divide). More info on perspective correction in the header of
% calcperspectiveCorrection.
% 
% 5. use fminsearch to find the matrix that minimizes the total distance
% between the target location and the apparent gaze location for each
% fixation.
% 
% References:
% http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.17.1329
% http://www.cse.psu.edu/~rtc12/CSE486/lecture12.pdf
% http://www.cse.psu.edu/~rtc12/CSE486/lecture13.pdf
% http://www.tomdalling.com/blog/modern-opengl/explaining-homogenous-coordinates-and-projective-geometry/


%% input parser

p = inputParser; p.KeepUnmatched = true;

% Required
p.addRequired('gazeDataFileName',@ischar);
p.addRequired('gazeCalFactorsFileName',@ischar);

% Optional analysis parameters
p.addParameter('fminsearchCalls',20,@isnumeric);
p.addParameter('testCalibration',true,@islogical);
p.addParameter('showFigures',false,@islogical);

% Optional display and I/O parameters
p.addParameter('verbosity','none', @ischar);

% Environment parameters
p.addParameter('tbSnapshot',[],@(x)(isempty(x) | isstruct(x)));
p.addParameter('timestamp',char(datetime('now')),@ischar);
p.addParameter('username',char(java.lang.System.getProperty('user.name')),@ischar);
p.addParameter('hostname',char(java.net.InetAddress.getLocalHost.getHostName),@ischar);

% parse
p.parse(gazeDataFileName, gazeCalFactorsFileName, varargin{:})


%% load gaze data
gazeData = load(gazeDataFileName);

targets.X = gazeData.gazeCalData.targets.X;
targets.Y = gazeData.gazeCalData.targets.Y;
pupil.X = gazeData.gazeCalData.pupil.X;
pupil.Y =gazeData.gazeCalData.pupil.Y;
glint.X = gazeData.gazeCalData.glint.X; 
glint.Y = gazeData.gazeCalData.glint.Y;
viewingDistance = gazeData.gazeCalData.viewingDistanceMm;

clear gazeData

%% Estimate the relative perspective correction factor

perspectiveCorrection = calcPerspectiveCorrection(targets,pupil,glint,viewingDistance);

%% Define apparent gaze vector and apply perspective correction (will be used in the fminsearch)
% theoretically, the Z term of the apparent gaze should include the
% distance between the pupil center and the center of the corneal curvature
% in pixels (is ranges between 50-200), and should look like this:
% 
%  apparentGaze = 1/perspectiveCorrection *[(Px -Gx); (Py -Gy); distPupilCornea - sqrt((Px -Gx)^2 +(Py -Gy)^2)];
% 
% Given the subsequent fminsearch step, we can set distancePupilCornea/perspectiveCorrection
% as 1 without loss of generality, so that:

%  apparentGaze = [(pX-gX)/perspectiveCorrection; (pY-gY)/perspectiveCorrection; (1 - sqrt(((pX-gX)/perspectiveCorrection)^2 + ((pY-gY)/perspectiveCorrection)^2))]

%% solve minimization problem
% initialize the matrix
X = [...
    1 0 0 0 ...
    0 1 0 0 ...
    0 0 0 0 ...
    0 0 0 1 ...
    ];

% exclude nan targets, if any
pupil.X = pupil.X(~isnan(targets.X));
pupil.Y = pupil.Y(~isnan(targets.X));
glint.X = glint.X(~isnan(targets.X));
glint.Y = glint.Y(~isnan(targets.X));
targets.X = targets.X(~isnan(targets.X));
targets.Y = targets.Y (~isnan(targets.X));

% Loop through calls to fminsearch, changing tolerance
for ff=1:p.Results.fminsearchCalls
    options = optimset('Display','off','MaxFunEvals', 10000,...
        'MaxIter', 10000, 'TolX',10^(-ff/2),'TolFun',10^(-ff/2),'PlotFcns',[] );
    [X, f] = fminsearch(@(param) ...
        errfun(param,pupil,glint,targets,viewingDistance,perspectiveCorrection),...
        X, options);
    if strcmp(p.Results.verbosity, 'full')
        disp(['RSS error: ',num2str(f)])
    end
end

% make the calibration matrix
calMat = [X(1:4); X(5:8); X(9:12); X(13:16)];


%% test the calibration if requested
if p.Results.testCalibration
    % get the calibrated fixations
    calibratedGaze.X = nan(size(pupil.X));
    calibratedGaze.Y = nan(size(pupil.X));
    tmp = nan(length(pupil.X),3);
    for ff = 1:length(pupil.X)
        pX = pupil.X(ff);
        pY = pupil.Y(ff);
        gX = glint.X(ff);
        gY = glint.Y(ff);
        aXYZW = calMat * [...
            (pX-gX)/perspectiveCorrection; ...
            (pY-gY)/perspectiveCorrection; ...
            (1 - sqrt(((pX-gX)/perspectiveCorrection)^2 + ((pY-gY)/perspectiveCorrection)^2)); ...
            1];
        tmp(ff,:) = (aXYZW(1:3)/aXYZW(4))';
    end
    calibratedGaze.X = tmp(:,1);
    calibratedGaze.Y = tmp(:,2);
    
    %get the calibration accuracy
    errors = sqrt((targets.X-calibratedGaze.X).^2+(targets.Y-calibratedGaze.Y).^2);
    accuracy = mean(errors(~isnan(errors)));
    
    if p.Results.showFigures
        % plot real target location and calibrated fixation location in screen
        % units and show the error
        figure; hold on;
        for ff = 1:length(pupil.X)
            plot(targets.X(ff), targets.Y(ff), 'bo')
            plot(calibratedGaze.X(ff), calibratedGaze.Y(ff),'rx')
            plot([targets.X(ff) calibratedGaze.X(ff)], [targets.Y(ff) calibratedGaze.Y(ff)],'k');
        end
        hold off;
        title(['Average error: ',num2str(accuracy),' mm. (viewing distance: ' num2str(viewingDistance) ' mm)'])
        xlabel('Horizontal position (mm)');ylabel('Vertical position (mm)');
        legend('Target position','Estimated gaze position')
    end
end


%% save out calibration matrix, perspectiveCorrection and metadata in a file
gazeCalFactors.calMatrix = calMat;
gazeCalFactors.perspectiveCorrection = perspectiveCorrection;
gazeCalFactors.meta = p.Results;
gazeCalFactors.meta.RSSerror = f;
if p.Results.testCalibration
    gazeCalFactors.meta.testedAccuracy = accuracy;
end
save(gazeCalFactorsFileName,'gazeCalFactors')

end % main function


%% error function for fminsearch
function errtot = errfun(param, pupil, glint, targets, viewingDistance, perspectiveCorrection)

err = nan(1,length(targets(:,1)));
CalMatrix = [param(1:4); param(5:8); param(9:12); param(13:16)];

% minimize error for each target
for i = 1:length(targets.X)
    pX = pupil.X(i);
    pY = pupil.Y(i);
    gX = glint.X(i);
    gY = glint.Y(i);
    x = targets.X(i);
    y = targets.Y(i);
    z = viewingDistance;
    
    % calibrated apparent gaze vector in homogeneous coordinates
    aXYZW = CalMatrix * [(pX-gX)/perspectiveCorrection; (pY-gY)/perspectiveCorrection; (1 - sqrt(((pX-gX)/perspectiveCorrection)^2 + ((pY-gY)/perspectiveCorrection)^2)); 1];
    
    % target vector in homogeneous coordinates
    oXYZW = [x; y; z; 1];
    
    % distance gaze - target in 3-D
    errXYZ = (aXYZW(1:3)/aXYZW(4)) - (oXYZW(1:3)/oXYZW(4));
    
    % error for this target
    err(i) = nansum(errXYZ.^2);
end

% total error
errtot = nansum(err.^2);
end % error function

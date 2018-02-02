function [calibratedPupil] = applySizeCalibration(pupilFileName,sizeCalFactorsFileName,varargin)
% Applies the size calibration factors to the pupil data.
%
% Description:
%   This routine applies the size calibration values to the pupil data.
% 
% Inputs:
%  pupilFileName              - name of the file with the pupil data to be calibrated, 
%                               as it results from the pupil pipeline.
%  sizeCalFactorsFileName     - name of the mat file to save the size
%       						conversion factor.
% 
% Optional key/value pairs (display and I/O)
%  'verbosity'                - level of verbosity. [none, full]
%
% Optional key/value pairs (analysis)
% 'calibratedPupilFileName'   - name of the output file containing the
%                               calibrated data, if the user wishes to save it on file.
%  'calibratedUnits'          - units in which the calibrated data is expressed
%                              (default [mm])
%  'whichFitToCalibrate'      - which of the pupil fit in the pupil file to calibrate
%   	                       (default 'radiusSmoothed', also available 'sceneConstrained' and 'initial')
%
% Optional key/value pairs (environment)
%  'tbSnapshot'               - the passed tbSnapshot output that is to be saved along
%      			                with the data
%  'timestamp'            	  - AUTOMATIC; The current time and date
%  'username'             	  - AUTOMATIC; The user
%  'hostname'             	  - AUTOMATIC; The host
%
% Outputs:
%  calibratedPupil            - structure containing the calibrated pupil width, height
%                               and area. The calibrated units are dependent on the size calibration
%                               method used.
% 

%% Parse vargin for options passed here
p = inputParser; p.KeepUnmatched = true;

% Required
p.addRequired('pupilFileName',@ischar);
p.addRequired('sizeCalFactorsFileName',@ischar);

% Optional analysis parameters
p.addParameter('calibratedPupilFileName','',@ischar);
p.addParameter('calibratedUnits','mm', @ischar);
p.addParameter('whichFitToCalibrate','radiusSmoothed', @ischar);

% Optional display and I/O parameters
p.addParameter('verbosity','none', @ischar);

% Environment parameters
p.addParameter('tbSnapshot',[],@(x)(isempty(x) | isstruct(x)));
p.addParameter('timestamp',char(datetime('now')),@ischar);
p.addParameter('username',char(java.lang.System.getProperty('user.name')),@ischar);
p.addParameter('hostname',char(java.net.InetAddress.getLocalHost.getHostName),@ischar);

% parse
p.parse(pupilFileName, sizeCalFactorsFileName, varargin{:})


%% load pupil data

tmpData = load(pupilFileName);
% pull transparent raw pupil data
rawPupilTransparent = tmpData.pupilData.(p.Results.whichFitToCalibrate).ellipses.values;
% convert to explicit data
for ii = 1 : size(rawPupilTransparent,1)
    rawPupilExplicit(ii,:) = ellipse_transparent2ex(rawPupilTransparent(ii,:));
end

% write pupilInPx as [Xaxis Yaxis area]
for ii = 1: size(rawPupilExplicit,1)
    % get the X and Y axis length according to orientation of the ellipse
    if round(cos(rawPupilExplicit(ii,5))) == 1
        horizontalAxis(ii) = rawPupilExplicit(ii,3) * 2;
        verticalAxis(ii) = rawPupilExplicit(ii,4) * 2;
    elseif round(cos(rawPupilExplicit(ii,5))) == 0
        horizontalAxis(ii) = rawPupilExplicit(ii,4) * 2;
        verticalAxis(ii) = rawPupilExplicit(ii,3) * 2;
    else
        horizontalAxis(ii) = NaN;
        verticalAxis(ii) = NaN;
    end
    % get ellipse area
    ellipseArea(ii) = rawPupilTransparent(ii,3);
end  % loop through frames
% gather the all raw values
pupilInPx = [horizontalAxis' verticalAxis' ellipseArea'];

clear tmpData
clear rawPupilTransparent
clear rawPupilExplicit


%% load calibration factors

tmpSizeCal = load(sizeCalFactorsFileName);
sizeCalFactors = tmpSizeCal.sizeCalFactors;
% check for warnings
if isfield(sizeCalFactors,'warnings')
    warning('There are some warnings for the size calibration factors, please check that the factors are legit');
end
% get the conversion factors
conversionFactors = [sizeCalFactors.PXperMM sizeCalFactors.PXperMM (sizeCalFactors.PXperMM).^2];

clear tmpSizeCal

%% apply conversion factors

calPupil = pupilInPx ./ conversionFactors;

%% save calibrated pupil and metadata
calibratedPupil.width = calPupil(:,1);
calibratedPupil.height = calPupil(:,2);
calibratedPupil.area = calPupil(:,3);
calibratedPupil.meta = p.Results;

if ~isempty(p.Results.calibratedPupilFileName)
    save(p.Results.calibratedPupilFileName, 'calibratedPupil');
end

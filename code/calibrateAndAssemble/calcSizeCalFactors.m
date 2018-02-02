function sizeCalFactors = calcSizeCalFactors(sizeDataFilesNames, sizeCalFactorsFileName, varargin)
% Calculates the factors needed for size calibration and scene
% geometry reconstruction.
%
%  >> DEV PLACEHOLDER: option to use scene geometry for size calibration
%  instead of size Videos will be implemented in this function
% 
% Description:
%   This routine computes the following calibration factors:
%     pxPerMm - linear size conversion factor derived from the ground truth
%       size of the circular calibration dots and the length of the major
%       axis of the fitted ellipses.
%     sceneDistanceMM{sceneDistancePX} - estimate of the scene distance in 
%       millimiters {pixels}, derived from the ground truth size of the
%       calibration dot, the size of the fitted ellipses and the camera
%       properties. If the camera properties are unknown, this will be
%       returned as empty.
%  
%   If more than a single size calibration dataset is used, the final
%   factors will be the mean of the factors obtained by each of datasets.
%   The routine automatically checks for the quality of the data
%   and returns warnings if the standard deviation exceedes a set of
%   arbitrary thresholds.
% 
% Inputs:
%  sizeDataFilesNames      - Cell array containing the names of the
%                            calibration data files to be used. Based on
%                            the selection of the size calibration type,
%                            this can be either a set of calibration dots
%                            video or a scene geometry file.
%  sizeFactorsFileName     - Name of the mat file to save the size
%                            conversion factor.
% 
% Optional key/value pairs:
%  'sizeCalibrationType'   - Type of size calibraton to be performed.
%  'sizeGroundTruthsInput' - Array containing the ground truth for the dot 
%                            size in mm. If left empty (default option),
%                            the routine will try to retrieve the ground
%                            truth from the sizeData files names following
%                            the instructions in the cell array
%                            groundTruthFinder.
%  'groundTruthFinder'     - Cell array with instruction to retrieve the 
%                            ground truths from the file name. It is
%                            composed as follows:
%                             { numOfDigits position referenceString} 
%                            where:
%                            numOfDigits = digits that compose the ground 
%                                          truth;
%                            position = either 'before' or 'after' 
%                                       the reference string.
%                            referenceString = key piece of string to 
%                                              locate the ground truth.
%                                              Note that it is case
%                                              insensitive.
%  'stdThreshold'           - 1x2 array with arbitrary thresholds for the 
%                             standard deviation of each calibration
%                             factor. If any of those values is exceeded,
%                             a warning is returned and saved with the
%                             output.
%  'cameraFocalLengthMM'    - Focal length of the camera in millimiters.
%  'cameraSensorSizeMM'     - Sensor size in mm in the format
%                             [HorizontalSizeMM VerticalSizeMM].
%  'sceneResolutionPX       - Resolution of the scene in pixels. Note that
%                             this might be different from the physical
%                             sensor resolution, and might be related to
%                             the compression/digitalization strategy
%                             applied to the video acquisition.
%  'cameraModel'            - Model to use to estimate the sceneDistance
%                             from the camera properties (default:
%                             'pinhole').
%
% Optional key/value pairs (display and I/O):
%  'verbosity'              - Level of verbosity. [none, full]
%
% Optional key/value pairs (environment)
%  'tbSnapshot'             - This should contain the output of the
%                             tbDeploymentSnapshot performed upon the
%                             result of the tbUse command. This documents
%                             the state of the system at the time of
%                             analysis.
%  'timestamp'              - AUTOMATIC; The current time and date
%  'username'               - AUTOMATIC; The user
%  'hostname'               - AUTOMATIC; The host
% 
% Outputs:
%   sizeCalFactors          - struct containing the conversion factors for pupil
%                             size in pixel per millimeters and the
%                             estimated scene distance in millimiters. In
%                             case the calibration is not accurate, a
%                             "warnings" field is created and saved to
%                             store warning texts about the accuracy
%                             problems. An additional meta field containing
%                             all input params is also saved.

%% Parse vargin for options passed here
p = inputParser; p.KeepUnmatched = true;

% Required
p.addRequired('sizeDataFilesNames',@(x) (iscell(x) | ischar(x)));
p.addRequired('sizeCalFactorsFileName',@ischar);

% Optional analysis parameters
p.addParameter('sizeCalibrationType','sizeVideos',@ischar) % alternative 'sceneGeometry'
p.addParameter('sizeGroundTruthsInput',[], @isnumeric)
p.addParameter('groundTruthFinder', {1 'before' 'mm'}, @iscell)
p.addParameter('stdThreshold', [0.5 8], @isnumeric)
p.addParameter('cameraFocalLengthMM', 16, @isnumeric)
p.addParameter('cameraSensorSizeMM', [4.69 3.54], @isnumeric)
p.addParameter('sceneResolutionPX', [640 480], @isnumeric)
p.addParameter('cameraModel', 'pinhole', @ischar)

% Optional display and I/O parameters
p.addParameter('verbosity','none', @ischar);

% Environment parameters
p.addParameter('tbSnapshot',[],@(x)(isempty(x) | isstruct(x)));
p.addParameter('timestamp',char(datetime('now')),@ischar);
p.addParameter('username',char(java.lang.System.getProperty('user.name')),@ischar);
p.addParameter('hostname',char(java.net.InetAddress.getLocalHost.getHostName),@ischar);

% parse
p.parse(sizeDataFilesNames, sizeCalFactorsFileName, varargin{:})

%% Check input and get ground truths

switch p.Results.sizeCalibrationType
    case 'sizeVideos'
% if the ground truths were input manually, check that they match the
% number of calibration data files.
if ~isempty (p.Results.sizeGroundTruthsInput)
    if iscell(sizeDataFilesNames) && length(sizeDataFilesNames)~=length(p.Results.sizeGroundTruthsInput)
        error('The number of data files and ground truth values does not match')
    elseif ~iscell(sizeDataFilesNames) && length(p.Results.sizeGroundTruthsInput) ~=1
        error('The number of data files and ground truth values does not match')
    else
        sizeGroundTruths = p.Results.sizeGroundTruthsInput;
    end
else
    % if no ground truths were input, try parsing the file name for ground
    % truths
    for rr = 1: length(sizeDataFilesNames)
        % get name only
        [~,RunName,~] = fileparts(sizeDataFilesNames{rr});
        % make it lower case
        lcRunName = lower(RunName);
        % make reference substring also lower case
        subString = lower(p.Results.groundTruthFinder{3});
        % find substring in the runName
        sIdx = strfind(lcRunName, subString);
        % switch position of the reference string
        switch p.Results.groundTruthFinder{2}
            case 'before'
                sizeGroundTruths(rr) = str2double(lcRunName(sIdx-p.Results.groundTruthFinder{1}:sIdx-1));
            case 'after'
                sizeGroundTruths(rr) = str2double(lcRunName(sIdx+length(subString):sIdx+length(subString)+p.Results.groundTruthFinder{1}));
        end
        % check that the retrieved sizeGroundTruths are legit
        if any(isnan(sizeGroundTruths))
            error('Some size ground truths could not be retrieved from the file names')
        end
    end
end


%% load all calibration data in a matrix
for rr = 1: length(sizeGroundTruths)
    % load in transparent form
    tmpData = load ([sizeDataFilesNames{rr} '_pupil.mat']);
    tmpTransparent = tmpData.pupilData.pPosteriorMeanTransparent;
    transparentData{rr} = tmpTransparent;
    % load in explicit form
    for ii = 1 : size(tmpTransparent,1)
        tmpExplicit(ii,:) = ellipse_transparent2ex(tmpTransparent(ii,:));
    end
    explicitData{rr} = tmpExplicit;
    clear tmpData
    clear tmpExplicit
    clear tmpTransparent
end

%% Derive conversion factors

% if camera properties are available, determine size of pixels in MM
if ~any([isempty(p.Results.cameraFocalLengthMM) isempty(p.Results.cameraSensorSizeMM) isempty(p.Results.sceneResolutionPX)])
    haveCameraFeatures = true;
    pixelSizeMM = p.Results.cameraSensorSizeMM./p.Results.sceneResolutionPX;
else
    haveCameraFeatures = false;
    pixelSizeMM = nan;
end

% loop through runs
for rr = 1: length(sizeGroundTruths)
    
    % loop through frames
    for ii = 1: size(explicitData{rr},1)
        
        % get the diameter of the calibration dot
        dotDiameterPX(ii) = explicitData{rr}(ii,3) * 2;
        
    end  % loop through frames
    
    % gather the all raw values
    rawValues{rr} = dotDiameterPX';
    
    % get an array for median value, std and conversion factor
    singleRunMedianPX(rr,:) = nanmedian(dotDiameterPX);
    singleRunStd(rr,:) = nanstd(dotDiameterPX);
    realSizeMM = sizeGroundTruths(rr);
    singleRunPXperMM(rr,:) = singleRunMedianPX(rr,:) ./ realSizeMM ;
    clear dotDiameterPX
    
    % if camera properties are available, use them to determine camera
    % distance (both in PX and MM)
    if haveCameraFeatures
        switch p.Results.cameraModel
            case 'pinhole'
                % convert median size for the single runs in MM
                singleRunMedianMM(rr,:) = singleRunMedianPX(rr,1)*pixelSizeMM(1);
                
                % apply the pinhole approximated formula to derive the
                % scene distance from the camera in MM
                singleRunSceneDistanceMM(rr,:) = realSizeMM * p.Results.cameraFocalLengthMM / singleRunMedianMM(rr,1);
                
                % also get the scene distance in pixels
                singleRunSceneDistancePX = singleRunSceneDistanceMM ./ pixelSizeMM(1);
                
            case 'realCamera'
                singleRunSceneDistanceMM(rr,:) = (p.Results.cameraFocalLengthMM * realSizeMM * p.Results.sceneResolutionPX(2))/(singleRunMedianPX(rr,1) *p.Results.cameraSensorSizeMM(2));
                
                % also get the scene distance in pixels
                singleRunSceneDistancePX = singleRunSceneDistanceMM ./ pixelSizeMM(1);
        end
    else
        singleRunSceneDistancePX = nan;
        singleRunSceneDistanceMM = nan;
    end
    
end % loop through runs


%% get the conversion factors as the mean of the individual ones
meanPXperMM = nanmean(singleRunPXperMM);
sdPXperMM = nanstd(singleRunPXperMM);

meanSceneDistancePX = nanmean(singleRunSceneDistancePX);
sdSceneDistancePX = nanstd(singleRunSceneDistancePX);

meanSceneDistanceMM = nanmean(singleRunSceneDistanceMM);
sdSceneDistanceMM = nanstd(singleRunSceneDistanceMM);

%% check if the calibration values are accurate
warningCounter = 0;
% check for standard deviation of the calibration factors
if  sdPXperMM > p.Results.stdThreshold(1)
    warningCounter = warningCounter +1;
    warningMessages{warningCounter} = 'High standard deviation for the calibration factors. The calibration might not be accurate.';
    warning(warningMessages{warningCounter})
end

    case 'sceneGeometry'
        error ('This option is not available yet')
end
%% compose sizeFactor struct

sizeCalFactors.PXperMM = meanPXperMM;
sizeCalFactors.sceneDistancePX = meanSceneDistancePX;
sizeCalFactors.sceneDistanceMM = meanSceneDistanceMM;

% add meta fields
sizeCalFactors.meta = p.Results;
sizeCalFactors.meta.sizeGroundTruths = sizeGroundTruths;
sizeCalFactors.meta.sdPXperMM = sdPXperMM;
sizeCalFactors.meta.sdSceneDistancePX = sdSceneDistancePX;
sizeCalFactors.meta.sdSceneDistanceMM = sdSceneDistanceMM;
if warningCounter > 0
    sizeCalFactors.warnings = warningMessages;
    clear warningMessage
end


%% save out data

save(sizeCalFactorsFileName,'sizeCalFactors')

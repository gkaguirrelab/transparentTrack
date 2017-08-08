function calcSizeFactors(sizeDataFilesNames, sizeFactorsFileName, varargin)
% calcSizeFactors(sizeDataFilesNames,sizeFactorsFileName)
%
% This routine computes the following size conversion factors:
%   horizontalPxPerMm
%   verticalPxPerMm
%   areaSqPxPerSqMm
% 
% If more than a single size calibration dataset is used, the size
% conversion factor will be the mean of the factors obtained by each of the
% size videos. The routine automatically checks for the quality of the data
% and returns warnings if the standard deviation exceedes a set of
% arbitrary thresholds.
%
% OUTPUTS:
%   sizeConversionFactor: 3 element vector of conversion factors for pupil
%   size expressed for the horizontal dim, vertical dim and area as
%   follows: [ hFactor vFactor aFactor]. The units are [ px/mm px/mm
%   sqPx/sqMm] respectively.
% 
% INPUTS:
%   dotDataFilesNames: cell array containing the names of the dot data
%       files to be used.
%   sizeFactorsFileName: name of the mat file to save the size
%       conversion factor.
% 
% Optional params:
%   sizeGroundTruths: array containing the ground truth for the dot size in
%       mm. If left empty (default option), the routine will try to
%       retrieve the ground truth from the sizeData files names following
%       the instructions in the cell array groundTruthFinder.
%   groundTruthFinder: cell array with instruction to retrieve the ground
%       truths from the file name. It is composed as follows:
%           { numOfDigits position referenceString} 
%           where:
%           numOfDigits = digits that compose the ground truth
%           position = either 'before' or 'after' the reference string
%           referenceString = key piece of string to locate the ground
%               truth. Note that it is case insensitive.
%   stdThreshold: 1x3 array with arbitrary thresholds for the standard
%       deviation of each calibration factor. If any of those values is
%       exceeded, a warning is returned and saved with the output.
%   areaErrorThreshold: arbitrary threshold value in squarePx per squareMm
%       to verify that the linear conversion factors are coherent with the
%       area conversion factor. If this threshold is exceeded a warning is
%       returned and saved with the output.
% 
% Optional key/value pairs (display and I/O)
%  'verbosity' - level of verbosity. [none, full]
%
% Options (environment)
%   tbSnapshot - the passed tbSnapshot output that is to be saved along
%      with the data
%   timestamp / username / hostname - these are automatically derived and
%      saved within the p.Results structure.
%
%% Parse vargin for options passed here
p = inputParser; p.KeepUnmatched = true;

% Required
p.addRequired('sizeDataFilesNames',@(x) (iscell(x) | ischar(x)));
p.addRequired('sizeFactorsFileName',@ischar);

% Optional analysis parameters
p.addParameter('sizeGroundTruths',[], @isnumeric)
p.addParameter('groundTruthFinder', {1 'before' 'mm'}, @iscell)
p.addParameter('stdThreshold', [0.5 0.5 4], @isnumeric)
p.addParameter('areaErrorThreshold', 2, @isnumeric)

% Optional display and I/O parameters
p.addParameter('verbosity','none', @ischar);

% Environment parameters
p.addParameter('tbSnapshot',[],@(x)(isempty(x) | isstruct(x)));
p.addParameter('timestamp',char(datetime('now')),@ischar);
p.addParameter('username',char(java.lang.System.getProperty('user.name')),@ischar);
p.addParameter('hostname',char(java.net.InetAddress.getLocalHost.getHostName),@ischar);

% parse
p.parse(sizeDataFilesNames, sizeFactorsFileName, varargin{:})

%% Check input and get ground truths

% if the ground truths were input manually, check that they match the
% number of calibration data files.
if ~isempty (p.Results.sizeGroundTruths)
    if iscell(sizeDataFilesNames) && length(sizeDataFilesNames)~=length(p.Results.sizeGroundTruths)
        error('The number of data files and ground truth values does not match')
    elseif ~iscell(sizeDataFilesNames) && length(p.Results.sizeGroundTruths) ~=1
        error('The number of data files and ground truth values does not match')
    else
        sizeGroundTruths = p.Results.sizeGroundTruths;
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
    tmpData = load (sizeDataFilesNames{rr});
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

% loop through runs
for rr = 1: length(sizeGroundTruths) 
    % loop through frames
   for ii = 1: size(explicitData{rr},1)
    % get the X and Y axis length according to orientation of the ellipse
    if round(cos(explicitData{rr}(ii,5))) == 1
        horizontalAxis(ii) = explicitData{rr}(ii,3) * 2;
        verticalAxis(ii) = explicitData{rr}(ii,4) * 2;
    elseif round(cos(explicitData{rr}(ii,5))) == 0
        horizontalAxis(ii) = explicitData{rr}(ii,4) * 2;
        verticalAxis(ii) = explicitData{rr}(ii,3) * 2;
    end
    % get ellipse area
    ellipseArea(ii) = transparentData{rr}(ii,3);
   end  % loop through frames
   % gather the all raw values
   rawValues{rr} = [horizontalAxis' verticalAxis' ellipseArea'];
   % get an array for median value, std and conversion factor
   singleRunMedian(rr,:) = [nanmedian(horizontalAxis) nanmedian(verticalAxis) nanmedian(ellipseArea)];
   singleRunStd(rr,:) = [nanstd(horizontalAxis) nanstd(verticalAxis) nanstd(ellipseArea)];
   realSize = [sizeGroundTruths(rr) sizeGroundTruths(rr) pi*((sizeGroundTruths(rr)/2)^2)];
   singleRunFactors(rr,:) = singleRunMedian(rr,:) ./ realSize ;
   clear horizontalAxis
   clear verticalAxis
   clear ellipseArea
end % loop through runs
    

%% get the conversion factors as the mean of the individual ones
sizeFactorsMean = mean(singleRunFactors);
sizeFactorsStd = std(singleRunFactors);


%% check if the calibration values are accurate
warningCounter = 0;
% check for standard deviation of the calibration factors
if any(sizeFactorsStd > p.Results.stdThreshold)
    warningCounter = warningCounter +1;
    warningMessages{warningCounter} = 'High standard deviation for the calibration factors. The calibration might not be accurate.';
    warning(warningMessages{warningCounter})
end

% check if linear factors and area factor are coherent
areaFromLinearFactors = sizeFactorsMean(1) * sizeFactorsMean(2);
areaFactor = sizeFactorsMean(3);
if abs(areaFactor - areaFromLinearFactors) > p.Results.areaErrorThreshold
    warningCounter = warningCounter +1;
    warningMessages{warningCounter} = 'The area conversion factor is not coherent with the linear conversion factors';
    warning(warningMessages{warningCounter})
end
%% compose sizeFactor struct

sizeFactors.horizontalPxPerMm = sizeFactorsMean(1);
sizeFactors.verticalPxPerMm = sizeFactorsMean(2);
sizeFactors.areaSqPxPerSqMm = sizeFactorsMean(3);

% add a meta field
sizeFactors.meta = p.Results;
if warningCounter > 0
    sizeFactors.warnings = warningMessages;
    clear warningMessage
end


%% save out data

save(sizeFactorsFileName,'sizeFactors')

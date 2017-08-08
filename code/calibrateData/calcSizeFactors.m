function calcSizeFactors(sizeDataFilesNames, sizeFactorsFileName, varargin)
% calcSizeFactors(sizeDataFilesNames,sizeFactorsFileName)
%
% This routine will compute the size conversion factor from px to mm.
% If more than a single size video is used, the size conversion factor will
% be the median of the factors obtained by each of the size videos.
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
%       mm. If left empty (default option), the routine will try to retrieve
%       the ground truth from the sizeData files names.

%
% 
% Optional key/value pairs (display and I/O)
%  'verbosity' - level of verbosity. [none, full]
%
%% Parse vargin for options passed here
p = inputParser; p.KeepUnmatched = true;

% Required
p.addRequired('dotDataFilesNames',@(x) (iscell(x) | ischar(x)));
p.addRequired('sizeConversionFileName',@ischar);

% Optional analysis parameters
p.addParameter('sizeGroundTruths',[], @isnumeric)
p.addParameter('groundTruthFinder', {1 'before' 'mm'}, @iscell)

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
        horizontalAxis(ii) = explicitData{rr}(ii,3);
        verticalAxis(ii) = explicitData{rr}(ii,4);
    elseif round(cos(explicitData{rr}(ii,5))) == 0
        horizontalAxis(ii) = explicitData{rr}(ii,4);
        verticalAxis(ii) = explicitData{rr}(ii,3);
    end
    % get ellipse area
    ellipseArea(ii) = transparentData{rr}(ii,3);
   end  % loop through frames
   % gather the all raw values
   sizeConversion.rawValues{rr} = [horizontalAxis' verticalAxis' ellipseArea'];
   % get an array for median value, std and conversion factor
   sizeConversion.singleRunMedian(rr,:) = [nanmedian(horizontalAxis) nanmedian(verticalAxis) nanmedian(ellipseArea)];
   sizeConversion.singleRunStd(rr,:) = [nanstd(horizontalAxis) nanstd(verticalAxis) nanstd(ellipseArea)];
   realSize = [sizeGroundTruths(rr) sizeGroundTruths(rr) pi*(sizeGroundTruths(rr)/2)^2];
   sizeConversion.singleRunFactors(rr,:) = sizeConversion.singleRunMedian(rr,:) ./ realSize ;
   clear horizontalAxis
   clear verticalAxis
   clear ellipseArea
end % loop through runs
    

%% get the conversion factors as the mean of the individual ones
sizeConversion.sizeFactorsMean = mean(sizeConversion.singleRunMedian);
sizeConversion.sizeFactorsStd = std(sizeConversion.sizeFactorsMean);


%% add a meta field
sizeConversion.meta = p.Results;


%% save out data
% add a meta field first
sizeConversion.meta = p.Results;
save(sizeFactorsFileName,'sizeConversion')

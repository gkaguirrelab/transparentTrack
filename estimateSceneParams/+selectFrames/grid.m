function [frameSet, gazeTargets, relativeCameraPosition] = grid(videoStemName, varargin)


%% input parser
p = inputParser; p.KeepUnmatched = true;

% Required
p.addRequired('videoStemName',@ischar);

p.addParameter('nBinsPerDimension',3,@isnumeric);
p.addParameter('badFrameErrorThreshold',2, @isnumeric);

% parse
p.parse(videoStemName, varargin{:})



% Load the timebase, pupilData, perimeter, and relative camera position for
% this acquisition
load([videoStemName '_pupil.mat'],'pupilData');
load([videoStemName '_relativeCameraPosition.mat'],'relativeCameraPosition');

% Extract the vectors to process
ellipses = pupilData.initial.ellipses.values;
ellipseFitRMSE = pupilData.initial.ellipses.RMSE;

% Identify the ellipses with RMSE fits that are below the "bad"
% threshold
goodFitIdx = find(ellipseFitRMSE < p.Results.badFrameErrorThreshold);
if isempty(goodFitIdx)
    error('No initial ellipse fits are good enough to guide the search; try adjusting badFrameErrorThreshold');
end

% First we divide the ellipse centers amongst a set of 2D bins across
% image space.
[ellipseCenterCounts,~,~,binXidx,binYidx] = ...
    histcounts2(ellipses(goodFitIdx,1),ellipses(goodFitIdx,2),p.Results.nBinsPerDimension);

% Anonymous functions for row and column identity given array position
rowIdx = @(b) fix( (b-1) ./ (size(ellipseCenterCounts,2)) ) +1;
colIdx = @(b) 1+mod(b-1,size(ellipseCenterCounts,2));

% Create a cell array of index positions corresponding to each of the
% 2D bins
idxByBinPosition = ...
    arrayfun(@(b) find( (binXidx==rowIdx(b)) .* (binYidx==colIdx(b)) ),1:1:numel(ellipseCenterCounts),'UniformOutput',false);

% Identify which bins are not empty
filledBinIdx = find(~cellfun(@isempty, idxByBinPosition));

% Identify the ellipse in each bin with the lowest fit RMSE
[~, idxMinErrorEllipseWithinBin] = arrayfun(@(x) nanmin(ellipseFitRMSE(goodFitIdx(idxByBinPosition{x}))), filledBinIdx, 'UniformOutput', false);
returnTheMin = @(binContents, x)  binContents(idxMinErrorEllipseWithinBin{x});
frameSet = cellfun(@(x) returnTheMin(goodFitIdx(idxByBinPosition{filledBinIdx(x)}),x),num2cell(1:1:length(filledBinIdx)));

% Create a set of nan gaze targets
gazeTargets = nan(2,length(frameSet));

% Return the relativeCameraPosition for these frames
relativeCameraPosition = relativeCameraPosition.values(:,frameSet);

end



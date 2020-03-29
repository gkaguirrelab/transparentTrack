function [frameSet, gazeTargets, relativeCameraPosition] = grid(videoStemName, varargin)


%% input parser
p = inputParser; p.KeepUnmatched = true;

% Required
p.addRequired('videoStemName',@ischar);

p.addParameter('nBinsPerDimension',5,@isnumeric);
p.addParameter('badFrameErrorThreshold',2, @isnumeric);
p.addParameter('minFramesPerBin',100, @isnumeric);

% parse
p.parse(videoStemName, varargin{:})



% Load the timebase, pupilData, perimeter, and relative camera position for
% this acquisition
load([videoStemName '_pupil.mat'],'pupilData');
load([videoStemName '_correctedPerimeter.mat'],'perimeter');
load([videoStemName '_relativeCameraPosition.mat'],'relativeCameraPosition');



%% Calculate likelhood SD across frames
% Obtain a measure for each frame of how completely the perimeter points
% define a full, 360 degrees around the pupil. This index of coverage of
% the pupil perimeter is distVals. If there is a perfectly uniform angular
% distribution of points in space around the pupil perimeter, then the
% distVals value will be zero. If there are perimeter points only at a
% single angular location around the pupil cirle, then the distVal will 1.

% The likelihood SD is based upon the RMSE of the fit of the elipse to the
% perimeter points for each frame
RMSE = pupilData.initial.ellipses.RMSE';

% Define the bins over which the distribution of perimeter angles will be
% evaluated. 20 bins works pretty well.
nDivisions = 20;
histBins = linspace(-pi,pi,nDivisions);

% Anonymous function returns the linear non-uniformity of a set of values,
% ranging from 0 when perfectly uniform to 1 when completely non-uniform.
nonUniformity = @(x) (sum(abs(x/sum(x)-mean(x/sum(x))))/2)/(1-1/length(x));

% Loop over frames. Frames which have no perimeter points will be given a
% distVal of NaN.
for ii = 1:length(RMSE)
    
    % Obtain the center of this fitted ellipse
    centerX = pupilData.initial.ellipses.values(ii,1);
    centerY = pupilData.initial.ellipses.values(ii,2);
    
    % Obtain the set of perimeter points
    Xp = perimeter.data{ii}.Xp;
    Yp = perimeter.data{ii}.Yp;
    
    % Calculate the deviation of the distribution of points from uniform
    linearNonUniformity(ii) = nonUniformity(histcounts(atan2(Yp-centerY,Xp-centerX),histBins));
end

% Subject the linearNonUniformity vector to a non-linear transformation.
% This has the effect of changing 0 --> 0.1, 0.8 --> 1, and values > 0.8
% --> infinity. This causes pupil perimeters with support at a single
% location (as opposed to fully around the perimeter) to have a markedly
% increased likelihood SD. Also, set InF values to something arbitrarily
% large.
distVals = (1./(1-sqrt(linearNonUniformity)))./10;
distVals(isinf(distVals)) = 1e20;

% The likelihood SD for each frame is the RMSE multiplied by the distVal
likelihoodPupilRadiusSDVector = distVals.*RMSE;

%% Find the best frame in each bin
ellipses = pupilData.initial.ellipses.values;
goodFitIdx = find(RMSE < p.Results.badFrameErrorThreshold);

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

% Identify the bins that have a sufficient number of frames to bother with
% looking for the best one
filledBinIdx = find(cellfun(@(x) size(x,1)>p.Results.minFramesPerBin, idxByBinPosition));

% Identify the ellipse in each bin with the lowest fit RMSE
[~, idxMinErrorEllipseWithinBin] = arrayfun(@(x) nanmin(likelihoodPupilRadiusSDVector(goodFitIdx(idxByBinPosition{x}))), filledBinIdx, 'UniformOutput', false);
returnTheMin = @(binContents, x)  binContents(idxMinErrorEllipseWithinBin{x});
frameSet = cellfun(@(x) returnTheMin(goodFitIdx(idxByBinPosition{filledBinIdx(x)}),x),num2cell(1:1:length(filledBinIdx)));

% Create a set of nan gaze targets
gazeTargets = nan(2,length(frameSet));

% Return the relativeCameraPosition for these frames
relativeCameraPosition = relativeCameraPosition.values(:,frameSet);


end



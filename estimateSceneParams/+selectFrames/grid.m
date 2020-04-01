function [frameSet, gazeTargets] = grid(videoStemName, varargin)
% Identify a set of frames that can be used to sync sceneGeometry
%
% Syntax:
%  [frameSet, gazeTargets] = grid(videoStemName)
%
% Description:
%   Positioning an eye model in a scene requires the selection of
%   informative frames of the acquisition to guide the alignment. This
%   routine selects frames that have a high-quality measurement of the
%   pupil perimeter, and are well distributed in gaze position.
%
% Inputs:
%	videoStemName         - Char vector. Full path to video file from which
%                           the scene observations have been derived. The
%                           stem name should omit the "_gray.avi" suffix
%                           that is usually present in the names of these
%                           video files.
%
% Optional key-value pairs:
%  'nBinsPerDimension'    - Scalar. Defines the number of divisions with
%                           which the ellipse centers are binned.
%  'badFrameErrorThreshold' - Scalar. Frames with RMSE values for the fit
%                           of an ellipse to the pupil perimeter above this
%                           threshold will not be selected to guide the
%                           scene parameter search.
%  'minFramesPerBin'      - Scalar. A given bin must have at least this
%                           many good frames for the best frame in the bin
%                           to be selected.
%
% Outputs:
%   frameSet              - A 1xm vector that specifies the m frame indices
%                           (indexed from 1) that identify the set of
%                           frames from the acquisition to guide the search
%   gazeTargets           - A 2xm matrix that provides the positions, in
%                           degrees of visual angle, of fixation targets
%                           that correspond to each of the frames. For this
%                           function, the returned matrix will be all nan
%                           values.
%

%% input parser
p = inputParser; p.KeepUnmatched = true;

% Required
p.addRequired('videoStemName',@ischar);

p.addParameter('nBinsPerDimension',8,@isnumeric);
p.addParameter('badFrameErrorThreshold',2, @isnumeric);
p.addParameter('minFramesPerBin',50, @isnumeric);

% parse
p.parse(videoStemName, varargin{:})


%% Load associated acquisition files

load([videoStemName '_pupil.mat'],'pupilData');
load([videoStemName '_correctedPerimeter.mat'],'perimeter');



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


end



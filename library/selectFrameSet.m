%  'nBinsPerDimension'    - Scalar. Defines the number of divisions with
%                           which the ellipse centers are binned.
%  'badFrameErrorThreshold' - Frames with RMSE fitting error above this
%                           threshold will not be selected to guide the
%                           scene parameter search.

if p.Results.verbose
        fprintf('Selecting ellipses to guide the search.\n');
    end
    
    % Identify the ellipses with RMSE fits that are below the "bad"
    % threshold
    goodFitIdx = find(ellipseFitRMSE < p.Results.badFrameErrorThreshold);
    if isempty(goodFitIdx)
        error('No initial ellipse fits are good enough to guide the search; try adjusting badFrameErrorThreshold');
    end
    
    % First we divide the ellipse centers amongst a set of 2D bins across
    % image space.
    [ellipseCenterCounts,Xedges,Yedges,binXidx,binYidx] = ...
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
    ellipseArrayList = cellfun(@(x) returnTheMin(goodFitIdx(idxByBinPosition{filledBinIdx(x)}),x),num2cell(1:1:length(filledBinIdx)));
    
    % If there is a fixationTargetArray, make sure it is empty
    if ~isempty(fixationTargetArray)
        warning('Cannot use fixationTargetArray unless ellipseArrayList is defined');
    end
    fixationTargetArray=[];
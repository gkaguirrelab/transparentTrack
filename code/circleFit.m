function [pCenters, pRadii,pMetric, gCenters, gRadii,gMetric, pupilRange, glintRange] = circleFit(I,params,pupilRange,glintRange)

% this function is used for both glint and pupil circle fitting.

%% set default thresholds for pupil and glint
% when using the function to find the glint the default pupil threshold
% will be used with no need to declare it, and viceversa.

if ~isfield(params,'pupilCircleThresh')
    params.pupilCircleThresh = 0.06;
end

if ~isfield(params,'glintCircleThresh')
    params.glintCircleThresh = 0.999;
end

if ~isfield(params,'glintOut')
    params.glintOut = 0.1;
end

%% circle fit

% create blurring filter
filtSize = round([0.01*min(params.imageSize) 0.01*min(params.imageSize) 0.01*min(params.imageSize)]);

% structuring element to dialate the glint
se = strel('disk',params.dilateGlint);

% Filter for pupil
padP = padarray(I,[size(I,1)/2 size(I,2)/2], 128);
h = fspecial('gaussian',[filtSize(1) filtSize(2)],filtSize(3));
pI = imfilter(padP,h);
pI = pI(size(I,1)/2+1:size(I,1)/2+size(I,1),size(I,2)/2+1:size(I,2)/2+size(I,2));
% Binarize pupil
binP = ones(size(pI));
binP(pI<quantile(double(pI(:)),params.pupilCircleThresh)) = 0;

% Filter for glint
gI = ones(size(I));
gI(I<quantile(double(pI(:)),params.glintCircleThresh)) = 0;
padG = padarray(gI,[size(I,1)/2 size(I,2)/2], 0);
h = fspecial('gaussian',[filtSize(1) filtSize(2)],filtSize(3));
gI = imfilter(padG,h);
gI = gI(size(I,1)/2+1:size(I,1)/2+size(I,1),size(I,2)/2+1:size(I,2)/2+size(I,2));
% Binarize glint
binG  = zeros(size(gI));
binG(gI>0.01) = 1;
dbinG = imdilate(binG,se);

% Silence the imfindcircles warning regarding circle size
warning('off','images:imfindcircles:warnForLargeRadiusRange');

% Find the pupil
[pCenters, pRadii,pMetric] = imfindcircles(binP,pupilRange,'ObjectPolarity','dark',...
    'Sensitivity',params.sensitivity);

% Find the glint
if ~params.pupilOnly
    [gCenters, gRadii,gMetric] = imfindcircles(dbinG,glintRange,'ObjectPolarity','bright',...
        'Sensitivity',params.sensitivity);
else
    gCenters = [NaN NaN];
    gRadii = NaN;
    gMetric = NaN;
end

% Restore the warning state
warning('on','images:imfindcircles:warnForLargeRadiusRange');

% Remove glints outside the pupil
if ~params.pupilOnly
    if ~isempty(pCenters) && ~isempty(gCenters)
        dists           = sqrt( (gCenters(:,1) - pCenters(1,1)).^2 + (gCenters(:,2) - pCenters(1,2)).^2 );
        gCenters(dists>(1 + params.glintOut)*(pRadii(1)),:) = [];
        gRadii(dists>(1 + params.glintOut)*(pRadii(1))) = [];
    end
end

% adjust the pupil range (for quicker processing)
if ~isempty(pCenters)
    pupilRange(1)   = min(floor(pRadii(1)*(1-params.rangeAdjust)),params.pupilRange(2));
    pupilRange(2)   = max(ceil(pRadii(1)*(1 + params.rangeAdjust)),params.pupilRange(1));
else
    pupilRange(1)   = max(ceil(pupilRange(1)*(1 - params.rangeAdjust)),params.pupilRange(1));
    pupilRange(2)   = min(ceil(pupilRange(2)*(1 + params.rangeAdjust)),params.pupilRange(2));
end

end
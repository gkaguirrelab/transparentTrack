function [pCenters, pRadii,pMetric, gCenters, gRadii,gMetric, pupilRange, glintRange] = circleFit(I,pupilCircleThresh,glintCircleThresh,pupilRange,glintRange,varargin)

% this function is used for both glint and pupil circle fitting.


%% parse input and define variables
p = inputParser;
% required input
p.addRequired('I');
p.addRequired('pupilCircleThresh', @isnumeric)
p.addRequired('glintCircleThresh', @isnumeric)
p.addRequired('pupilRange',@isnumeric);
p.addRequired('glintRange',@isnumeric);


% optional inputs
pupilOnlyDefault = false;
glintOutDefault = 0.1;
dilateGlintDefault = 6;
sensitivityDefault = 0.99;
rangeAdjustDefault = 0.05;
p.addParameter('pupilOnly', pupilOnlyDefault, @islogical);
p.addParameter('glintOut', glintOutDefault, @isnumeric);
p.addParameter('dilateGlint', dilateGlintDefault, @isnumeric);
p.addParameter('sensitivity', sensitivityDefault, @isnumeric);
p.addParameter('rangeAdjust', rangeAdjustDefault, @isnumeric);

%parse
p.parse(I,pupilCircleThresh,glintCircleThresh,pupilRange,glintRange,varargin{:})

% define optional variables values
pupilOnly = p.Results.pupilOnly;
glintOut = p.Results.glintOut;
dilateGlint = p.Results.dilateGlint;
sensitivity = p.Results.sensitivity;
rangeAdjust = p.Results.rangeAdjust;

%% circle fit

% create blurring filter
filtSize = round([0.01*min(size(I)) 0.01*min(size(I)) 0.01*min(size(I))]);

% structuring element to dialate the glint
se = strel('disk',dilateGlint);

% Filter for pupil
padP = padarray(I,[size(I,1)/2 size(I,2)/2], 128);
h = fspecial('gaussian',[filtSize(1) filtSize(2)],filtSize(3));
pI = imfilter(padP,h);
pI = pI(size(I,1)/2+1:size(I,1)/2+size(I,1),size(I,2)/2+1:size(I,2)/2+size(I,2));
% Binarize pupil
binP = ones(size(pI));
binP(pI<quantile(double(pI(:)),pupilCircleThresh)) = 0;

% Filter for glint
gI = ones(size(I));
gI(I<quantile(double(pI(:)),glintCircleThresh)) = 0;
padG = padarray(gI,[size(I,1)/2 size(I,2)/2], 0);
h = fspecial('gaussian',[filtSize(1) filtSize(2)],filtSize(3));
gI = imfilter(padG,h);
gI = gI(size(I,1)/2+1:size(I,1)/2+size(I,1),size(I,2)/2+1:size(I,2)/2+size(I,2));
% Binarize glint
binG  = zeros(size(gI));
binG(gI>0.01) = 1;
dbinG = imdilate(binG,se);

% store the warning state
origWarnState = warning;

% Silence the imfindcircles warning regarding circle size
warning('off','images:imfindcircles:warnForLargeRadiusRange');
warning('off','images:imfindcircles:warnForSmallRadius');

% Find the pupil
[pCenters, pRadii,pMetric] = imfindcircles(binP,pupilRange,'ObjectPolarity','dark',...
    'Sensitivity',sensitivity);

% Find the glint
if ~pupilOnly
    [gCenters, gRadii,gMetric] = imfindcircles(dbinG,glintRange,'ObjectPolarity','bright',...
        'Sensitivity',sensitivity);
else
    gCenters = [NaN NaN];
    gRadii = NaN;
    gMetric = NaN;
end

% Restore the warning state
warning(origWarnState);

% Remove glints outside the pupil
if ~pupilOnly
    if ~isempty(pCenters) && ~isempty(gCenters)
        dists = sqrt( (gCenters(:,1) - pCenters(1,1)).^2 + (gCenters(:,2) - pCenters(1,2)).^2 );
        gCenters(dists>(1 + glintOut)*(pRadii(1)),:) = [];
        gRadii(dists>(1 + glintOut)*(pRadii(1))) = [];
    end
end

% adjust the pupil range (for quicker processing)
if ~isempty(pCenters)
    pupilRange(1)   = min(floor(pRadii(1)*(1-rangeAdjust)),pupilRange(2)); %%% CHECH THIS, it was params.pupilRange 
    pupilRange(2)   = max(ceil(pRadii(1)*(1 + rangeAdjust)),pupilRange(1));
else
    pupilRange(1)   = max(ceil(pupilRange(1)*(1 - rangeAdjust)),pupilRange(1));
    pupilRange(2)   = min(ceil(pupilRange(2)*(1 + rangeAdjust)),pupilRange(2));
end

end
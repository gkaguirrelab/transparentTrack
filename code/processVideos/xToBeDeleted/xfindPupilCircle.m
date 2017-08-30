function [pCenters, pRadii,pMetric, pupilRange] = findPupilCircle(I,pupilCircleThresh,pupilRange,imfindcirclesSensitivity,rangeAdjust)
% findGlintAndPupilCircles(I,pupilCircleThresh,glintCircleThresh,pupilRange,glintRange,pupilOnly,glintOut,dilateGlint,imfindcirclesSensitivity,rangeAdjust)
%
% this function is used for pupil circle fitting.


%% parse input and define variables
p = inputParser;
% required input
p.addRequired('I');
p.addRequired('pupilCircleThresh', @isnumeric)
p.addRequired('pupilRange',@isnumeric);
p.addRequired('imfindcirclesSensitivity', @isnumeric);
p.addRequired('rangeAdjust', @isnumeric);

% parse
p.parse(I,pupilCircleThresh,pupilRange,imfindcirclesSensitivity,rangeAdjust);

%% circle fit

% create blurring filter
filtSize = round([0.01*min(size(I)) 0.01*min(size(I)) 0.01*min(size(I))]);


% Filter for pupil
padP = padarray(I,[size(I,1)/2 size(I,2)/2], 128);
h = fspecial('gaussian',[filtSize(1) filtSize(2)],filtSize(3));
pI = imfilter(padP,h);
pI = pI(size(I,1)/2+1:size(I,1)/2+size(I,1),size(I,2)/2+1:size(I,2)/2+size(I,2));
% Binarize pupil
binP = ones(size(pI));
binP(pI<quantile(double(pI(:)),pupilCircleThresh)) = 0;


% store the warning state
origWarnState = warning;

% Silence the imfindcircles warning regarding circle size
warning('off','images:imfindcircles:warnForLargeRadiusRange');
warning('off','images:imfindcircles:warnForSmallRadius');

% Find the pupil
[pCenters, pRadii,pMetric] = imfindcircles(binP,pupilRange,'ObjectPolarity','dark',...
    'Sensitivity',imfindcirclesSensitivity);

% Restore the warning state
warning(origWarnState);


% adjust the pupil range (for quicker processing)
if ~isempty(pCenters)
    pupilRange(1)   = min(floor(pRadii(1)*(1-rangeAdjust)),pupilRange(2));  
    pupilRange(2)   = max(ceil(pRadii(1)*(1 + rangeAdjust)),pupilRange(1));
else
    pupilRange(1)   = max(ceil(pupilRange(1)*(1 - rangeAdjust)),pupilRange(1));
    pupilRange(2)   = min(ceil(pupilRange(2)*(1 + rangeAdjust)),pupilRange(2));
end

end % main function
function [lidFitParams] = findEyelidBounds(thisFrame, irisEllipseFitTransparent )

%% Parse input and define variables
p = inputParser;

% required input
p.addRequired('thisFrame',@isnumeric);
p.addRequired('irisEllipseFitTransparent',@isnumeric);

% Optional analysis params
p.addParameter('ringIdxThreshold', 2, @isnumeric);
p.addParameter('irisRingExpandFraction', 0.2, @isnumeric);
p.addParameter('dipThresholdWithinBand', 0.04, @isnumeric);
p.addParameter('dipThresholdBetweenBands', 0.02, @isnumeric);
p.addParameter('pX0UpperLidTransparent', [160, 260, 100000, .8, pi],@isnumeric);
p.addParameter('pX0LowerLidTransparent', [160, 60, 100000, .8, pi],@isnumeric);
p.addParameter('UpperLidTransparentUB',[160, 300, 100000, .8, pi],@isnumeric);
p.addParameter('UpperLidTransparentLB',[160, 100, 100000, .8, pi],@isnumeric);
p.addParameter('LowerLidTransparentUB',[160, 150, 100000, .8, pi],@isnumeric);
p.addParameter('LowerLidTransparentLB',[160, 0, 100000, .8, pi],@isnumeric);



%% Parse and check the parameters
p.parse(thisFrame, irisEllipseFitTransparent)

% Obtain the iris and sclera bands, based upon dilating and contracting the
% elliptical fit to the iris border
videoSizeX = size(thisFrame,1);
videoSizeY = size(thisFrame,2);
[Ye, Xe] = ind2sub(size(thisFrame),1:1:videoSizeX*videoSizeY);
expandedIrisEllipseTransparent = irisEllipseFitTransparent';
contractedIrisEllipseTransparent = irisEllipseFitTransparent';
expandedIrisEllipseTransparent(3) = expandedIrisEllipseTransparent(3)*(1+p.Results.irisRingExpandFraction);
contractedIrisEllipseTransparent(3) = contractedIrisEllipseTransparent(3)*(1-p.Results.irisRingExpandFraction);
[~,tmpEllipseEdgeDistance] = ellipse_distance(Xe, Ye, ellipse_ex2im(ellipse_transparent2ex(expandedIrisEllipseTransparent)));
[scleraBandY, scleraBandX] = ind2sub(size(thisFrame), find(abs(tmpEllipseEdgeDistance) < p.Results.ringIdxThreshold ));
[~,tmpEllipseEdgeDistance] = ellipse_distance(Xe, Ye, ellipse_ex2im(ellipse_transparent2ex(contractedIrisEllipseTransparent)));
[irisBandY, irisBandX] = ind2sub(size(thisFrame), find(abs(tmpEllipseEdgeDistance) < p.Results.ringIdxThreshold ));

% Define an array of scalars that adjust the parameters of the ellipse
% definition so that a minimum change of 1 unit is acceptable for each param
paramScaler = [1 1 0.001 100 10];

% Define the objective function
myObjFun = @(x) hartigansWeightedBandCount(x, ...
    irisBandX, irisBandY, ...
    scleraBandX, scleraBandY, ...
    thisFrame, ...
    paramScaler);

% Combine the upper and lower lid parameters into single vectors
x0 = [p.Results.pX0UpperLidTransparent p.Results.pX0LowerLidTransparent];
ub = [p.Results.UpperLidTransparentUB p.Results.LowerLidTransparentUB];
lb = [p.Results.UpperLidTransparentLB p.Results.LowerLidTransparentLB];

% Apply the param
x0 = x0 .* [paramScaler paramScaler];
ub = ub .* [paramScaler paramScaler];
lb = lb .* [paramScaler paramScaler];

% Set some options
options = optimset('fmincon');
options = optimset(options,'Diagnostics','off','Display','iter-detailed', ...
    'LargeScale','off','Algorithm','sqp', 'DiffMinChange',1);

% Fit that sucker
[lidFitParams,e] = fmincon(myObjFun, x0, [], [], [], [], lb, ub, [], options);

% Scale the area parameter back up
lidFitParams = lidFitParams ./ [paramScaler paramScaler];

% Display the fit
imshow(thisFrame)
hold on

pFitImplicit = ellipse_ex2im(ellipse_transparent2ex(lidFitParams(1:5)'));
fh=@(x,y) pFitImplicit(1).*x.^2 +pFitImplicit(2).*x.*y +pFitImplicit(3).*y.^2 +pFitImplicit(4).*x +pFitImplicit(5).*y +pFitImplicit(6);
fimplicit(fh,[1, max([videoSizeX videoSizeY]), 1, max([videoSizeX videoSizeY])],'Color', 'green','LineWidth',1.5);

pFitImplicit = ellipse_ex2im(ellipse_transparent2ex(lidFitParams(6:10)'));
fh=@(x,y) pFitImplicit(1).*x.^2 +pFitImplicit(2).*x.*y +pFitImplicit(3).*y.^2 +pFitImplicit(4).*x +pFitImplicit(5).*y +pFitImplicit(6);
fimplicit(fh,[1, max([videoSizeX videoSizeY]), 1, max([videoSizeX videoSizeY])],'Color', 'green','LineWidth',1.5);

drawnow


end % findEyelidBounds function



function [objVal] = hartigansWeightedBandCount(lidEllipseParams, irisBandX, irisBandY, scleraBandX, scleraBandY, thisFrame, paramScaler)
% function [ceq] = evalHartigansCriteria(lidEllipseParams, irisBandX, irisBandY, scleraBandX, scleraBandY, thisFrame, dipThresholdWithinBand, dipThresholdBetweenBands, paramScaler)
%
% This function implements a non-linear constraint that requires that
%  1) the Hartigan's dip test is less than the dipThresholdWithinBand for
%     both the irisBandValues and the scleraBandValues. This is an
%     indication that the bands contain only these image components, and
%     not eyelid.
%  2) the Hartigan's dip test is greater than the dipThresholdBetweenBands
%     for a vector that combines the iris and sclera band valyes. This is
%     an indication that there are separate pools of dark and light pixels.

% Scale the area ellipse parameters back up
lidEllipseParams = lidEllipseParams ./ [paramScaler paramScaler];

% unpack the lidEllipseParams into the upper and lower lids. Each lid is
% defined by the edge of an ellipse cast in transparent form
upperLidEllipseParams = lidEllipseParams(1:5);
lowerLidEllipseParams = lidEllipseParams(6:10);

% identify the subset of iris points that are within the
% region of the palpebral fissure, and then obtain the values within the
% image frame for thes epoints
[~,tmpEllipseEdgeDistanceUpper] = ellipse_distance(irisBandX, irisBandY, ellipse_ex2im(ellipse_transparent2ex(upperLidEllipseParams)));
[~,tmpEllipseEdgeDistanceLower] = ellipse_distance(irisBandX, irisBandY, ellipse_ex2im(ellipse_transparent2ex(lowerLidEllipseParams)));
inBoundIrisPoints = find((tmpEllipseEdgeDistanceUpper < 0) .* (tmpEllipseEdgeDistanceLower < 0));
irisBandValues = double(thisFrame(sub2ind(size(thisFrame),irisBandY(inBoundIrisPoints),irisBandX(inBoundIrisPoints))));

% repeat this for the sclera values
[~,tmpEllipseEdgeDistanceUpper] = ellipse_distance(scleraBandX, scleraBandY, ellipse_ex2im(ellipse_transparent2ex(upperLidEllipseParams)));
[~,tmpEllipseEdgeDistanceLower] = ellipse_distance(scleraBandX, scleraBandY, ellipse_ex2im(ellipse_transparent2ex(lowerLidEllipseParams)));
inBoundScleraPoints = find((tmpEllipseEdgeDistanceUpper < 0) .* (tmpEllipseEdgeDistanceLower < 0));
scleraBandValues = double(thisFrame(sub2ind(size(thisFrame),scleraBandY(inBoundScleraPoints),scleraBandX(inBoundScleraPoints))));

% get the log of the number of pixels in the bands
numPointsScore = log10(length(inBoundIrisPoints) + length(inBoundScleraPoints));

% obtain the Hartigan's Dip criteria for the iris and scleral band values;
% we want these values to be small, reflecting a unimodal set of values
dipTestIrisBand = HartigansDipTest(sort(irisBandValues));
dipTestScleraBand = HartigansDipTest(sort(scleraBandValues));

objVal = -1 * (numPointsScore / (10^(dipTestIrisBand*100)) / (10^(dipTestScleraBand*100)));

end % evalHartigansCriteria



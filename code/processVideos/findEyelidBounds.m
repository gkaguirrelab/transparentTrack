function [  ] = findEyelidBounds(thisFrame, irisEllipseFitTransparent )

%% Parse input and define variables
p = inputParser;

% required input
p.addRequired('thisFrame',@isnumeric);
p.addRequired('irisEllipseFitTransparent',@isnumeric);

% Optional analysis params
p.addParameter('ringIdxThreshold', 2, @isnumeric);
p.addParameter('irisRingExpandFraction', 0.2, @isnumeric);
p.addParameter('pX0UpperLidTransparent', [160, 250, 100000, .8, pi],@isnumeric);
p.addParameter('pX0LowerLidTransparent', [160, 60, 100000, .8, pi],@isnumeric);

p.addParameter('ellipseTransparentLB',[0, 0, 10000, 0, -0.5*pi],@isnumeric);
p.addParameter('ellipseTransparentUB',[240,320,30000,0.2, 0.5*pi],@isnumeric);


%% Parse and check the parameters
p.parse(thisFrame, irisEllipseFitTransparent)

videoSizeX = size(thisFrame,1);
videoSizeY = size(thisFrame,2);
[Ye, Xe] = ind2sub(size(thisFrame),1:1:videoSizeX*videoSizeY);
expandedIrisEllipseTransparent = irisEllipseFitTransparent';
contractedIrisEllipseTransparent = irisEllipseFitTransparent';
expandedIrisEllipseTransparent(3) = expandedIrisEllipseTransparent(3)*(1+p.Results.irisRingExpandFraction);
contractedIrisEllipseTransparent(3) = contractedIrisEllipseTransparent(3)*(1-p.Results.irisRingExpandFraction);
[~,tmpEllipseEdgeDistance] = ellipse_distance(Xe, Ye, ellipse_ex2im(ellipse_transparent2ex(expandedIrisEllipseTransparent)));
[Ysclera, Xsclera] = ind2sub(size(thisFrame), find(abs(tmpEllipseEdgeDistance) < p.Results.ringIdxThreshold ));
[~,tmpEllipseEdgeDistance] = ellipse_distance(Xe, Ye, ellipse_ex2im(ellipse_transparent2ex(contractedIrisEllipseTransparent)));
[Yiris, Xiris] = ind2sub(size(thisFrame), find(abs(tmpEllipseEdgeDistance) < p.Results.ringIdxThreshold ));

[~,tmpEllipseEdgeDistanceUpper] = ellipse_distance(Xiris, Yiris, ellipse_ex2im(ellipse_transparent2ex(p.Results.pX0UpperLidTransparent)));
[~,tmpEllipseEdgeDistanceLower] = ellipse_distance(Xiris, Yiris, ellipse_ex2im(ellipse_transparent2ex(p.Results.pX0LowerLidTransparent)));

[Yiris, Xiris] = ind2sub(size(thisFrame), find( (tmpEllipseEdgeDistanceUpper < 0) .* (tmpEllipseEdgeDistanceLower < 0) ));


[dip] = HartigansDipTest(sort(thisFrame(sub2ind(size(thisFrame),Yiris,Xiris))));




end


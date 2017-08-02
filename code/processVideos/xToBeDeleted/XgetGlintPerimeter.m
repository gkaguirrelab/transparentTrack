function [binG] = getGlintPerimeter (I, gCenters, gRadii, glintEllipseThresh)

% create a mask from circle fitting parameters (note: glint
% is already dilated
glintMask = zeros(size(I));
glintMask = insertShape(glintMask,'FilledCircle',[gCenters(1,1) gCenters(1,2) gRadii(1)],'Color','white');
glintMask = im2bw(glintMask);

% apply mask to grey image
maskedGlint = immultiply(I,glintMask);

% convert back to gray
gI = uint8(maskedGlint);

% Binarize glint
binG  = ones(size(gI));
binG(gI<quantile(double(I(:)),glintEllipseThresh)) = 0;

% get perimeter of glint
binG = bwperim(binG);
function [binP] = getPupilPerimeter(I,pCenters,pRadii, sep, params)
% get binary pupil perimeter. This function will apply a pre-determined
% circular mask to the frame, binarize the resulting "patch" image with a
% user determined threshold, and extract the perimeter of the bigger region
% surviving the thresholding process (believed to be the pupil).

% create a mask from circle fitting parameters
pupilMask = zeros(size(I));
pupilMask = insertShape(pupilMask,'FilledCircle',[pCenters(1,1) pCenters(1,2) pRadii(1)],'Color','white');
pupilMask = imdilate(pupilMask,sep);
pupilMask = im2bw(pupilMask);

% apply mask to grey image complement image
cI = imcomplement(I);
maskedPupil = immultiply(cI,pupilMask);

% convert back to gray
pI = uint8(maskedPupil);
% Binarize pupil
binP = ones(size(pI));
binP(pI<quantile(double(cI(:)),params.ellipseThresh(1))) = 0;

% remove small objects
binP = bwareaopen(binP, 500);

% fill the holes
binP = imfill(binP,'holes');

% get perimeter of object
binP = bwperim(binP);
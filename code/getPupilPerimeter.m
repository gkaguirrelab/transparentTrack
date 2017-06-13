function [binP] = getPupilPerimeter(I,pCenters,pRadii, pupilEllipseThresh, varargin)
% get binary pupil perimeter. This function will apply a pre-determined
% circular mask to the frame, binarize the resulting "patch" image with a
% user determined threshold, and extract the perimeter of the bigger region
% surviving the thresholding process (believed to be the pupil).

%% parse input and define variables

p = inputParser;
% required input
p.addRequired('I');
p.addRequired('pCenters');
p.addRequired('pRadii');
p.addRequired('pupilEllipseThresh');

% optional inputs
maskBoxDefault = [4 30];
p.addParameter('maskBox', maskBoxDefault, @isnumeric);

%parse
p.parse(I,pCenters,pRadii, pupilEllipseThresh, varargin{:})

% define optional variables values
maskBox = p.Results.maskBox;

%% create a mask from circle fitting parameters

% structuring element for pupil mask size
sep = strel('rectangle',maskBox);
        
% generate mask
pupilMask = zeros(size(I));
pupilMask = insertShape(pupilMask,'FilledCircle',[pCenters(1,1) pCenters(1,2) pRadii(1)],'Color','white');
pupilMask = imdilate(pupilMask,sep);
pupilMask = im2bw(pupilMask);

%% apply mask to grey image complement image
cI = imcomplement(I);
maskedPupil = immultiply(cI,pupilMask);

%% convert back to gray
pI = uint8(maskedPupil);

%% Binarize pupil
binP = ones(size(pI));
binP(pI<quantile(double(cI(:)),pupilEllipseThresh)) = 0;

%% remove small objects
binP = bwareaopen(binP, 500);

%% fill the holes
binP = imfill(binP,'holes');

%% get perimeter of object
binP = bwperim(binP);
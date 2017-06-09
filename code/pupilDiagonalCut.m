function binPcut = pupilDiagonalCut (I, cutoutU, cutoutR, binP, underGlint, overGlint, leftGlint, rightGlint)
% this function cuts the pupil perimeter diagonally. To be used with
% trackPupil.

[Xp, Yp] = ind2sub(size(binP),find(binP));

[~, idx] = sort(rightGlint);

% get the cut perimeter (upper cut)
binPcutU = zeros(size(I));
binPcutU(sub2ind(size(binP),Xp(underGlint),Yp(underGlint))) = 1;
binPcutU(sub2ind(size(binP),Xp(overGlint),Yp(overGlint))) = 1;
binPcutU(sub2ind(size(binP),Xp(overGlint(idx(1+round(cutoutU/2):end - round(cutoutU/2)))),Yp(overGlint(idx(1+round(cutoutU/2):end - round(cutoutU/2)))))) = 0;

% get the cut perimeter (right  cut)
binPcutR = zeros(size(I));
binPcutR(sub2ind(size(binP),Xp(leftGlint),Yp(leftGlint))) = 1;
binPcutR(sub2ind(size(binP),Xp(rightGlint),Yp(rightGlint))) = 1;
binPcutR(sub2ind(size(binP),Xp(rightGlint(idx(cutoutR:end))),Yp(rightGlint(idx(cutoutR:end))))) = 0;

% remove small objects
binPcutM = immultiply(binPcutU,binPcutR);
CC = bwconncomp(binPcutM);
numPixels = cellfun(@numel,CC.PixelIdxList);
[~,idx] = max(numPixels);
binPcut = zeros(size(I));
binPcut(CC.PixelIdxList{idx}) = 1;
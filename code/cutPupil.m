function [binPcut] = cutPupil (binP,U,R,Xg,Yg)

% cuts the pupil perimeter according to the instructions in U and R, where:
% U = distance in pixels from the glint on the Y axis, going UP
% R = distance in pixels from the glint on the X axis going RIGHT
%
% If R is INF, the routine will perform an horizontal cut, thus removing
% the upper part of the pupil perimeter, U pixels above the glint.
% If U is INF, the routine will perform a vertical cut, thus removing the
% right part of the pupil, R pixel from the glint.
% If R and U are both finite integers, a combinations of the 2 cuts will be
% performed, effectively performing a diagonal cut of only the pixels
% common to the 2 orthogonal cuts.
%
%  The function accepts negative values for U and R, making the cut go
%  below/left of the glint, but it is not recommended to cut more than half
%  of the perimeter.

%% get perimeter points
[Yp, Xp] = ind2sub(size(binP),find(binP));

%% store full perimeter length
fullPerimeterLength = length(Xp);

%% find pixels

% only horizontal cut
if isinf(R) && ~isinf(U)
    underGlint = find (Yp > Yg - U);  %%% MUST CHECK IF EMPTY
    % get the cut perimeter
    binPcut = zeros(size(binP));
    binPcut(sub2ind(size(binP),Yp(underGlint),Xp(underGlint))) = 1;
    % remove small objects
    CC = bwconncomp(binPcut);
    numPixels = cellfun(@numel,CC.PixelIdxList);
    [biggest,idx] = max(numPixels);
    if idx >0
        binPcut = zeros(size(binP));
        binPcut(CC.PixelIdxList{idx}) = 1;
    end
    
% only vertical cut
elseif isinf(U) && ~isinf(R)
    leftOfGlint = find (Xp < Xg + R);  %%% MUST CHECK IF EMPTY
    % get the cut perimeter
    binPcut = zeros(size(binP));
    binPcut(sub2ind(size(binP),Yp(leftOfGlint),Xp(leftOfGlint))) = 1;
    % remove small objects
    CC = bwconncomp(binPcut);
    numPixels = cellfun(@numel,CC.PixelIdxList);
    [biggest,idx] = max(numPixels);
    if idx >0
        binPcut = zeros(size(binP));
        binPcut(CC.PixelIdxList{idx}) = 1;
    end
    
% diagonal cut
elseif ~isinf(U) && ~isinf(R)
    underGlint = find (Yp > Yg - U);
    leftOfGlint = find (Xp < Xg + R);
    toKeep = intersect(underGlint, leftOfGlint);
    % get the cut perimeter
    binPcut = zeros(size(binP));
    binPcut(sub2ind(size(binP),Yp(toKeep),Xp(toKeep))) = 1;
    % remove small objects
    CC = bwconncomp(binPcut);
    numPixels = cellfun(@numel,CC.PixelIdxList);
    [biggest,idx] = max(numPixels);
    if idx >0
        binPcut = zeros(size(binP));
        binPcut(CC.PixelIdxList{idx}) = 1;
    end
    
% return an error
else
    error('U and R can''t be INF at the same time')

end

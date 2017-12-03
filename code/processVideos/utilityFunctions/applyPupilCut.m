function [binPcut] = applyPupilCut(binP, radiusThresh, theta)
% Trim the pupil perimeter in an image by a specified amount
%
% Description:
%   Removes points from the boundary of pupil according to the passed
%   parameters. First, the non-zero points in the passed image are found.
%   The center of these points is taken as the reference point. Then,
%   points that are more distant from the center than a threshold line are
%   removed. Finally, the points are returned to the image.
%
% Inputs:
%   binP          - The image that contains the pupil perimeter
%   radiusThresh  - The distance (in pixels) from the reference point that
%                   defines the extent of boundary point removal
%   theta         - The angle (radians) with respect to the reference point
%                   of the threshold line. Convention:
%                       0     - left
%                       pi/2  - up
%                       pi    - right
%                       3pi/2 - down
%
% Output:
%   binPcut       - Image containing the refined pupil boundary
%


% initialize output matrix
binPcut = zeros(size(binP));

% create an anonymous function to return a rotation matrix given theta in
% radians
returnRotMat = @(theta) [cos(theta) -sin(theta); sin(theta) cos(theta)];

% get the perimeter points
[Yp, Xp] = ind2sub(size(binP),find(binP));

% create a matrix that will center the pupil boundary points relative to
% the center of the boundary points
xCenter=min(Xp)+(max(Xp)-min(Xp))/2;
yCenter=min(Yp)+(max(Yp)-min(Yp))/2;
centerMatrix = repmat([xCenter'; yCenter'], 1, length(Xp));

% Rotate the boundary points by theta radians
rotatedPoints = ...
    returnRotMat(theta) * ([Xp,Yp]' - centerMatrix) + centerMatrix;

% identify those points that have a position that exceeds radiusThresh in
% the Y axis of the rotated points
toBeKeptIdx = (find(rotatedPoints(2,:) <= (yCenter+radiusThresh)))';

% place the retained points in the output matrix
if ~isempty(toBeKeptIdx)
    linearInd = sub2ind(size(binP), Yp(toBeKeptIdx), Xp(toBeKeptIdx));
    binPcut(linearInd)=255;
end

end % function
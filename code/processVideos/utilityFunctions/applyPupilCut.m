function [binPcut] = applyPupilCut(binP,radiusThresh,theta)
% applyPupilCut(binP,radiusThresh,theta)
%
% cuts the pupil perimeter according to the instructions in radius and
% theta where:
%   radius = distance in pixels from the center of the perimeter
%   theta = polar angle, ranging from zero (superior vertical) to 2*pi
%

binPcut = binP*0;

% Create an anonymous function to return a rotation matrix given theta in
% radians
returnRotMat = @(theta) [cos(theta) -sin(theta); sin(theta) cos(theta)];

%% get perimeter points
[Yp, Xp] = ind2sub(size(binP),find(binP));

% Create a matrix that will center the pupil boundary points relative to
% the center of the boundary points
xCenter=min(Xp)+(max(Xp)-min(Xp))/2;
yCenter=min(Yp)+(max(Yp)-min(Yp))/2;
centerMatrix = repmat([xCenter'; yCenter'], 1, length(Xp));

% Rotate the boundary points by theta radians
rotatedPoints = returnRotMat(theta) * ([Xp,Yp]' - centerMatrix) + centerMatrix;
                
% Identify those points that have a position that exceeds radiusThresh in
% the Y axis of the rotated points
toBeKeptIdx = (find(rotatedPoints(2,:) <= (yCenter+radiusThresh)))';

% Image the kept points
if ~isempty(toBeKeptIdx)
    linearInd = sub2ind(size(binP), Yp(toBeKeptIdx), Xp(toBeKeptIdx));
    binPcut(linearInd)=255;
end


end % function
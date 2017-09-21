function perspectiveCorrection = calcPerspectiveCorrection(targets,pupil,glint,viewingDistance)
% perspectiveCorrection = calcPerspectiveCorrection(targets,pupil,glint,viewingDistance)
% 
% This function calculates a plausible term to operate homogeneous divide
% on the apparent gaze vector, so that it is "projected" on the same plane
% as the targets (the screen viewed at a known viewing distance).
% 
% Let's assume that the target plane (screen) and the apparent gaze plane
% share the same center. Given that the apparent gaze is defined as
% (pupil-glint) and that the light source that generates the glint and the
% camera hole are almost coincident, we can say that when apparentGaze = 0,
% the subject is looking directly into the camera. 
% 
% We then pick a "target far from the center of the screen" (roughly the
% center of the gaze plane) to make sure that the corresponding apparent
% gaze vector is non-zero (and likely among the longest available for this
% calibration).
% 
% We then have a target vector defined as:
% 
%     targetVector = [Tx Ty viewingDistance] 
% 
% and an apparent gaze vector defined as:
% 
%     apparentGazeVector = [(Px -Gx) (Py - Gy) Z], with Z unknown.
% 
% By simple geometric proportion:
% 
%  (Tx^2 + Tx^2 + viewingDist^2)/(Tx^2 + Tx^2) = ((Px -Gx)^2 + (Py - Gy)^2 + Z^2) /((Px -Gx)^2 + (Py - Gy)^2)
%  
% Then:
%   
%  sqrt(((Px -Gx)^2 + (Py - Gy)^2 + Z^2)) = sqrt((Tx^2 + Tx^2 + viewingDist^2)/(Tx^2 + Tx^2)) * sqrt((Px -Gx)^2 + (Py - Gy)^2)
% 
% Without loss in generality we can assume Z >> (Px -Gx) 
%  and Z >> (Py -Gy)^2, and simplify the above formula in:
%  
%  Z = sqrt((Tx^2 + Tx^2 + viewingDist^2)/(Tx^2 + Tx^2)) * sqrt((Px -Gx)^2 + (Py - Gy)^2)
% 
% This is our estimate for the relative perspective correction factor.
% 
% 
% MORE ON PERSPECTIVE CORRECTION
% 
% To project the apparent gaze on the same plane as the target we need to
% operate the following perspective correction matrix:
% 
%   [1 0 0 0 ...
%    0 1 0 0 ...
%    0 0 1 0 ...                            
%    0 0 1 0 ]     * [(Px -Gx); (Py - Gy); Z; 1]
% 
%  and the following homogeneous divide
%       [ perspectiveCorrectedGaze ] = 1/Z * [(Px -Gx); (Py - Gy); Z; Z]
% 


%% Pull out variables
targets = [targets.X, targets.Y];
pupil = [pupil.X, pupil.Y];
glint = [glint.X, glint.Y];
%% Find the corner target furthest from the center
[~,I] = max(nansum(abs(targets),2));
cT = targets(I, :);
cP = pupil(I, :);
cG = glint(I, :);
%% Calculate the perspective correction factor
perspectiveCorrection = (sqrt((cT(1))^2 + (cT(2))^2 + viewingDistance^2) ...
    / sqrt((cT(1))^2 + (cT(2))^2)) * sqrt((cP(1) - cG(1))^2 + (cP(2) - cG(2))^2);
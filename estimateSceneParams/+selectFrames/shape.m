function [frameSet, gazeTargets] = shape(videoStemName, rhoTarget, thetaTarget)
% Identify a fixation frame that can be used to sync sceneGeometry
%
% Syntax:
%  [frameSet, gazeTargets] = selectFrames.shape(videoStemName)
%
% Description:
%   Positioning an eye model in a scene requires the selection of
%   informative frames of the acquisition to guide the alignment.
%
%   On some occasions, we lack a measurement of gaze pose. This routine is
%   finds a video frame that contains an image of the pupil that is most
%   similar in shape (independent of size) to a target shape.
%
% Inputs:
%	videoStemName         - Char vector. Full path to video file from which
%                           the scene observations have been derived. The
%                           stem name should omit the "_gray.avi" suffix
%                           that is usually present in the names of these
%                           video files.
%
% Optional key-value pairs:
%   none
%
% Outputs:
%   frameSet              - Scalar that specifies a frame index
%                           (indexed from 1).
%   gazeTargets           - A 2x1 matrix that provides the positions, in
%                           degrees of visual angle of fixation positions
%                           for the frames. The "shape" selection process
%                           returns "nan" values for this variable.
%


%% input parser
p = inputParser; p.KeepUnmatched = true;

% Required
p.addRequired('videoStemName',@ischar);
p.addRequired('rhoTarget',@isscalar);
p.addRequired('thetaTarget',@isscalar);
p.addParameter('distValsThreshold',0.275, @isnumeric);

% parse
p.parse(videoStemName, rhoTarget, thetaTarget)


%% Load files
load([videoStemName '_pupil.mat'],'pupilData');
load([videoStemName '_correctedPerimeter.mat'],'perimeter');
load([videoStemName '_glint.mat'],'glintData');


%% Obtain the theta and rho values
rho = pupilData.initial.ellipses.values(:,4);
rho = 1-sqrt(1-rho.^2);
theta = pupilData.initial.ellipses.values(:,5);
theta = theta.*2;


%% Calculate the quality of match to the target theta and rho

% This vector expresses the difference in pupil shape between the reference
% period of the sceneGeometryIn and the frames of the acquisition.
shapeDifference = ...
    sqrt(rhoTarget^2 + rho.^2 - 2*rhoTarget.*rho.*cos(thetaTarget-theta))./2;

% This error vector will be minimized for frames on which the shape of the
% pupil is most similar to the reference period of the sceneGeometryIn.
shapeMatchError = sqrt(sum(shapeDifference.^2,2));


%% Calculate likelhood SD across frames
% Obtain a measure for each frame of how completely the perimeter points
% define a full, 360 degrees around the pupil. This index of coverage of
% the pupil perimeter is distVals. If there is a perfectly uniform angular
% distribution of points in space around the pupil perimeter, then the
% distVals value will be zero. If there are perimeter points only at a
% single angular location around the pupil cirle, then the distVal will 1.

% The likelihood SD is based upon the RMSE of the fit of the elipse to the
% perimeter points for each frame
RMSE = pupilData.initial.ellipses.RMSE';

% Define the bins over which the distribution of perimeter angles will be
% evaluated. 20 bins works pretty well.
nDivisions = 20;
histBins = linspace(-pi,pi,nDivisions);

% Anonymous function returns the linear non-uniformity of a set of values,
% ranging from 0 when perfectly uniform to 1 when completely non-uniform.
nonUniformity = @(x) (sum(abs(x/sum(x)-mean(x/sum(x))))/2)/(1-1/length(x));

% Loop over frames. Frames which have no perimeter points will be given a
% distVal of NaN.
for ii = 1:length(RMSE)
    
    % Obtain the center of this fitted ellipse
    centerX = pupilData.initial.ellipses.values(ii,1);
    centerY = pupilData.initial.ellipses.values(ii,2);
    
    % Obtain the set of perimeter points
    Xp = perimeter.data{ii}.Xp;
    Yp = perimeter.data{ii}.Yp;
    
    % Calculate the deviation of the distribution of points from uniform
    linearNonUniformity(ii) = nonUniformity(histcounts(atan2(Yp-centerY,Xp-centerX),histBins));
end

% Subject the linearNonUniformity vector to a non-linear transformation.
% This has the effect of changing 0 --> 0.1, 0.8 --> 1, and values > 0.8
% --> infinity. This causes pupil perimeters with support at a single
% location (as opposed to fully around the perimeter) to have a markedly
% increased likelihood SD. Also, set InF values to something arbitrarily
% large.
distVals = (1./(1-sqrt(linearNonUniformity)))./10;
distVals(isinf(distVals)) = 1e20;
distVals(isnan(distVals)) = 1e20;

% Adopt a threshold above which a partial pupil perimeter will not be used
distVals(distVals>p.Results.distValsThreshold) = 1e20;

% Set values with no glint to an arbitrarily large number
distVals(isnan(glintData.X)) = 1e20;

% The likelihood SD for each frame is the RMSE multiplied by the distVal
likelihoodPupilRadiusSDVector = (distVals.*RMSE)';


%% Find the best frame
% The best frame will have an acceptable likelihoodPupilRadiusSDVector, and
% the lowest shapeMatchError
[~, bestFrame] = nanmin(likelihoodPupilRadiusSDVector .* shapeMatchError);

% Return the best frame
frameSet = bestFrame;

% We don't really know the gazeTarget for this frame. The calling routine
% might know, but we will set the value for now to nan.
gazeTargets = [nan; nan];

end



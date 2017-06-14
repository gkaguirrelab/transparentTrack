function [glint, glintTrackingParams] = trackGlint(grayI, glintFile, varargin)

% This function tracks the glint using the circle patch + direct ellipse
% fitting approach.
% 
% There usually is no need to change the parameters for glint tracking, as
% it is pretty consistently tracked with the default settings.
% 
% 
% Output
% ======
%       glint file, glint variable, glintTrackingParams.
% 
% Input params
% ============
%       grayI : 3D array of gray frames to track
%       glintFile : name of the output matFile in which to save the glint
%         results.
%       
% Options
% =======
%       displayTracking : display online glint tracking in a figure
%           (default: false)
%       gammaCorrection : gamma correction to be applied in current frame
%       glintCircleThresh : threshold value to locate the glint for circle
%           fitting (default 0.999)
%       glintRange : radius range for cirfle fitting of the glint (default [10 30])
%       glintEllipseThresh : threshold value to locate the glint for
%           ellipse fitting (default 0.9)
% 
% 
% Usage example
% =============
%  [glint, glintTrackingParams] = trackGlint(grayI, glintFile, 'displayTracking', true)



%% parse input and define variables

p = inputParser;
% required input
p.addRequired('grayI');
p.addRequired('glintFile',@isstr);

% optional inputs
displayTrackingDefault = false;
gammaCorrectionDefault = 1;
glintCircleThreshDefault =  0.999;
glintRangeDefault = [10 30];
glintEllipseThreshDefault = 0.9;
p.addParameter('displayTracking', displayTrackingDefault, @islogical);
p.addParameter('gammaCorrection', gammaCorrectionDefault, @isnumeric);
p.addParameter('glintCircleThresh', glintCircleThreshDefault, @isnumeric);
p.addParameter('glintRange', glintRangeDefault, @isnumeric);
p.addParameter('glintEllipseThresh', glintEllipseThreshDefault, @isnumeric);

%parse
p.parse(grayI, glintFile, varargin{:})

% define optional variables values
displayTracking = p.Results.displayTracking;
gammaCorrection = p.Results.gammaCorrection;
glintCircleThresh =  p.Results.glintCircleThresh;
glintRange = p.Results.glintRange;
glintEllipseThresh = p.Results.glintEllipseThresh;


%% Initialize glint struct

% display main tracking parameters
disp('Starting tracking with the following parameters:');
disp('Circle threshold: ')
disp(glintCircleThresh)
disp('Ellipse threshold: ')
disp(glintEllipseThresh)

% get number of frames from grayI
numFrames = size(grayI,3);

% main glint params
glint.X = nan(numFrames,1);
glint.Y = nan(numFrames,1);
glint.size = nan(numFrames,1);

% glint fit params
glint.implicitEllipseParams = nan(numFrames,6);
glint.explicitEllipseParams= nan(numFrames,5);
glint.distanceErrorMetric= nan(numFrames,1);

% glint mask params
glint.circleRad = nan(numFrames,1);
glint.circleX = nan(numFrames,1);
glint.circleY = nan(numFrames,1);
glint.circleStrength = nan(numFrames,1);

% glint flags
glint.ellipseFittingError = nan(numFrames,1);

%% Track the glint

% NOTE: it is necessary to give a starting value for the pupil range,
% to remove the glint outside the pupil during the circle fitting step.
% The value will be updated during the tracking and pupil results won't
% be stored at this point.
pupilRange = [30 90];
pupilCircleThresh = 0.06;

if displayTracking
    ih = figure;
end

%loop through frames
for i = 1:numFrames 
    % Get the frame
    I = squeeze(grayI(:,:,i));
    
    % adjust gamma for this frame
    I = imadjust(I,[],[],gammaCorrection);
    
    % Show the frame (optional)
    if displayTracking
            imshow(I, 'Border', 'tight')
    end
    
    % track with circles (using the default options)
    [~,~,~, gCenters, gRadii,gMetric, pupilRange, glintRange] = circleFit(I,pupilCircleThresh,glintCircleThresh,pupilRange,glintRange);
    
    % get a more precise tracking with direct ellipse fitting
    if isempty(gCenters) %no glint was found by circleFit
        continue
    else % glint was found by circleFit
        % getGlintPerimeter
        [binG] = getGlintPerimeter (I, gCenters, gRadii, glintEllipseThresh);
        % Fit ellipse to glint
        [Xg, Yg] = ind2sub(size(binG),find(binG));
        try
            Egi = ellipsefit_direct(Xg,Yg);
            Eg = ellipse_im2ex(Egi);
            % get errorMetric
            [~,dg,~,~] = ellipse_distance(Xg, Yg, Egi);
            gdistanceErrorMetric = nanmedian(sqrt(sum(dg.^2)));
        catch ME
        end
        if  exist ('ME', 'var')
            glint.X(i)= gCenters(1,1);
            glint.Y(i) = gCenters(1,2);
            glint.size(i) = gRadii(1);
            glint.circleStrength(i) = gMetric(1);
            glint.ellipseFittingError(i) = 1;
            clear ME
        end
        
        % store results
        if exist ('Eg','var')
            if ~isempty (Eg) && isreal(Egi)
                glint.X(i) = Eg(2);
                glint.Y(i) = Eg(1);
                glint.circleStrength(i) = gMetric(1);
                glint.implicitEllipseParams(i,:) = Egi';
                glint.explicitEllipseParams(i,:) = Eg';
                glint.distanceErrorMetric(i) = gdistanceErrorMetric;
                % circle params for glint
                glint.circleStrength(i) = gMetric(1);
                glint.circleRad(i) = gRadii(1);
                glint.circleX(i) = gCenters(1,1);
                glint.circleY(i) = gCenters(1,2);
            end
            clear Eg Egi errors
        else
            glint.X(i)= gCenters(1,1);
            glint.Y(i) = gCenters(1,2);
            glint.size(i) = gRadii(1);
            glint.circleStrength(i) = gMetric(1);
        end
    end
    
    % plot results
    if displayTracking && ~isnan(glint.X(i))
        hold on
        plot(glint.X(i),glint.Y(i),'+b');
        hold off
    end
end

close all

% store the tracking params
glintTrackingParams = p.Results;

%% save out a mat file with the glint tracking data
save (glintFile, 'glint', 'glintTrackingParams')
    

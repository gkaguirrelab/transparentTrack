function [glint, params] = trackGlint(grayI, glintFile, varargin)

% This function tracks the glint using the circle patch + direct ellipse
% fitting approach.
% 
% There usually is no need to change the parameters for glint tracking, as
% it is pretty consistently tracked with the default settings.
% 
% Input params
%       grayI = series of gray frames to track
%       glintFile = name of the matFile in which to save the glint
%         results.
%       
% Options (can be modified if needed)
%       displayTracking : display online glint tracking in a figure
%           (default: false)
%       glintCircleThresh
%       glintRange
%       glintOut
%       glintEllipseThresh
%       gammaCorrection
% 
% Additional options (it is advised not to modify them)
%       rangeAdjust
%       sensitivity
% 
% 
% 
% Output: glint file, glint variable, params for control file.



%% parse input and define variables

p = inputParser;
% required input
p.addRequired('grayI',isa(grayI,'uint8'));
p.addRequired('glintFile',@isstr);

% optional inputs
displayTrackingDefault = false;
gammaCorrectionDefault = 1;
glintCircleThreshDefault =  0.999;
glintRangeDefault = [10 30];
glintEllipseThreshDefault = 0.9;
p.addParameter('displayTracking', displayTrackingDefault, @islogic);
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




% %% set default params
% 
% % params for circle patch
% if ~isfield(params,'rangeAdjust')
%     params.rangeAdjust = 0.05;
% end
% if ~isfield(params,'circleThresh')
%     params.glintCircleThresh =  0.999;
% end
% if ~isfield(params,'glintRange')
%     params.glintRange = [10 30];
% end
% if ~isfield(params,'glintOut')
%     params.glintOut = 0.1;
% end
% if ~isfield(params,'sensitivity')
%     params.sensitivity = 0.99;
% end
% if ~isfield(params,'dilateGlint')
%     params.dilateGlint = 5;
% end
% 
% % params for direct ellipse fitting
% if ~isfield(params,'ellipseThresh')
%     params.glintEllipseThresh = 0.9;
% end
% 
% 
% % param to display online tracking
% if ~isfield(params,'displayTracking')
%     params.displayTracking = 0;
% end

%% Initialize glint struct

% display main tracking parameters
disp('Starting tracking with the following parameters:');
disp('Circle threshold: ')
disp(glintCircleThresh)
disp('Ellipse threshold: ')
disp(params.glintEllipseThresh)

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
pupilCircleThresh = 0.6;

%loop through frames
for i = 1:numFrames 
    % Get the frame
    I = squeeze(grayI(:,:,i));
    
    % adjust gamma for this frame
    I = imadjust(I,[],[],gammaCorrection);
    
    % Show the frame (optional)
    if displayTracking
        imshow(I);
    end
    
    % track with circles (using the default options)
    [~,~,~, gCenters, gRadii,gMetric, pupilRange, glintRange] = circleFit(I,pupilCircleThresh,glintCircleThresh,pupilRange,glintRange);
    
    % get a more precise tracking with direct ellipse fitting
    if isempty(gCenters) %no glint was found by circleFit
        continue
    else % glint was found by circleFit
        % getGlintPerimeter
        [binG] = getGlintPerimeter (I, gCenters, gRadii, params);
        
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
        else
            glint.X(i)= gCenters(1,1);
            glint.Y(i) = gCenters(1,2);
            glint.size(i) = gRadii(1);
            glint.circleStrength(i) = gMetric(1);
        end
        if params.displayTracking && ~isnan(glint.X(i))
            plot(glint.X(i),glint.Y(i),'+b');
        end
        clear Eg Egi errors
    end
end
close I

%% save out a mat file with the glint tracking data
save (glintFile, 'glint', 'params')
    

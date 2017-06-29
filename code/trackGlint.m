function [glint, glintTrackingParams] = trackGlint(grayVideoName, glintFileName, varargin)

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
%       grayVideoName : name and path of the gray video on which to track
%           the glint
%       glintFileName : name of the output matFile in which to save the glint
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
%  [glint, glintTrackingParams] = trackGlint(grayVideoName, glintFileName, 'displayTracking', true)



%% parse input and define variables

p = inputParser;
% required input
p.addRequired('grayVideoName',@isstr);
p.addRequired('glintFile',@isstr);

% optional control parameters
p.addParameter('displayTracking', false, @islogical);
p.addParameter('gammaCorrection', 1, @isnumeric);
p.addParameter('glintCircleThresh', 0.999, @isnumeric);
p.addParameter('glintRange', [10 30], @isnumeric);
p.addParameter('glintEllipseThresh', 0.9, @isnumeric);

% Environment parameters
p.addParameter('tbSnapshot',[],@(x)(isempty(x) | isstruct(x)));
p.addParameter('timestamp',char(datetime('now')),@ischar);
p.addParameter('hostname',char(java.lang.System.getProperty('user.name')),@ischar);
p.addParameter('username',char(java.net.InetAddress.getLocalHost.getHostName),@ischar); 

%parse
p.parse(grayVideoName, glintFileName, varargin{:})

% define optional variables values
displayTracking = p.Results.displayTracking;
gammaCorrection = p.Results.gammaCorrection;
glintCircleThresh =  p.Results.glintCircleThresh;
glintRange = p.Results.glintRange;
glintEllipseThresh = p.Results.glintEllipseThresh;



%% read video file
% load pupilPerimeter
inObj = VideoReader(grayVideoName);

% get number of frames
numFrames = floor(inObj.Duration*inObj.FrameRate);



%% Initialize glint struct

% display main tracking parameters
disp('Starting tracking with the following parameters:');
disp('Circle threshold: ')
disp(glintCircleThresh)
disp('Ellipse threshold: ')
disp(glintEllipseThresh)


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


% initialize progress bar
progBar = ProgressBar(numFrames,'Tracking the glint...');

%loop through frames
for ii = 1:numFrames 
    % Get the frame
    I = readFrame(inObj);
    
    % adjust gamma for this frame
    I = imadjust(I,[],[],gammaCorrection);
    
    I = rgb2gray (I);
    % track with circles (using the default options)
    [~,~,~, gCenters, gRadii,gMetric, pupilRange, glintRange] = circleFit(I,pupilCircleThresh,glintCircleThresh,pupilRange,glintRange);
    
    % get a more precise tracking with direct ellipse fitting
    if isempty(gCenters) %no glint was found by circleFit
        if ~mod(ii,10);progBar(ii);end % update progressbar
        continue
    else % glint was found by circleFit
        % getGlintPerimeter
        [binG] = getGlintPerimeter (I, gCenters, gRadii, glintEllipseThresh);
        % Fit ellipse to glint
        [Xg, Yg] = ind2sub(size(binG),find(binG));
        try
            % turn of warnings for singular matrix
            origWarnState = warning();
            warning('off','MATLAB:singularMatrix');
            warning('off','MATLAB:illConditionedMatrix');
            Egi = ellipsefit_direct(Xg,Yg);
            warning(origWarnState);
            Eg = ellipse_im2ex(Egi);
            % get errorMetric
            [~,dg,~,~] = ellipse_distance(Xg, Yg, Egi);
            gdistanceErrorMetric = nanmedian(sqrt(sum(dg.^2)));
        catch ME
        end
        if  exist ('ME', 'var')
            glint.X(ii)= gCenters(1,1);
            glint.Y(ii) = gCenters(1,2);
            glint.size(ii) = gRadii(1);
            glint.circleStrength(ii) = gMetric(1);
            glint.ellipseFittingError(ii) = 1;
            clear ME
        end
        
        % store results
        if exist ('Eg','var')
            if ~isempty (Eg) && isreal(Egi)
                glint.X(ii) = Eg(2);
                glint.Y(ii) = Eg(1);
                glint.circleStrength(ii) = gMetric(1);
                glint.implicitEllipseParams(ii,:) = Egi';
                glint.explicitEllipseParams(ii,:) = Eg';
                glint.distanceErrorMetric(ii) = gdistanceErrorMetric;
                % circle params for glint
                glint.circleStrength(ii) = gMetric(1);
                glint.circleRad(ii) = gRadii(1);
                glint.circleX(ii) = gCenters(1,1);
                glint.circleY(ii) = gCenters(1,2);
            end
            clear Eg Egi errors
        else
            glint.X(ii)= gCenters(1,1);
            glint.Y(ii) = gCenters(1,2);
            glint.size(ii) = gRadii(1);
            glint.circleStrength(ii) = gMetric(1);
        end
        if ~mod(ii,10);progBar(ii);end % update progressbar
    end
end

close all
clear inObj

% store the tracking params
glintTrackingParams = p.Results;

%% save out a mat file with the glint tracking data
save (glintFileName, 'glint','glintTrackingParams')
    

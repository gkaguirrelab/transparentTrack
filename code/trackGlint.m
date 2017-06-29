function [glintData] = trackGlint(grayVideoName, glintFileName, varargin)
% function [glintData] = trackGlint(grayVideoName, glintFileName, varargin)
%
%   DO WE NEED TO SAVE ALL THIS GLINT TRACKING DATA?
%       Most of it seems inessential. Also, there is row x column
%       disagrement in the output.
%
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
%       gammaCorrection : gamma correction to be applied in current frame
%       glintCircleThresh : threshold value to locate the glint for circle
%           fitting (default 0.999)
%       glintRange : radius range for cirfle fitting of the glint (default [10 30])
%       glintEllipseThresh : threshold value to locate the glint for
%           ellipse fitting (default 0.9)
%       pupilRange: DESCRIBE WHAT THIS IS
%       pupilCircleThresh: DESCRIBE WHAT THIS IS
%
% Usage example
% =============
%  [glint, glintTrackingParams] = trackGlint(grayVideoName, glintFileName, 'displayTracking', true)



%% parse input and define variables

p = inputParser;
% required input
p.addRequired('grayVideoName',@isstr);
p.addRequired('glintFileName',@isstr);

% optional control parameters
p.addParameter('gammaCorrection', 1, @isnumeric);
p.addParameter('glintCircleThresh', 0.999, @isnumeric);
p.addParameter('glintRange', [10 30], @isnumeric);
p.addParameter('glintEllipseThresh', 0.9, @isnumeric);
p.addParameter('pupilRange', [30 90], @isnumeric);
p.addParameter('pupilCircleThresh', 0.06, @isnumeric);
p.addParameter('verbosity', 'none', @ischar);

% Environment parameters
p.addParameter('tbSnapshot',[],@(x)(isempty(x) | isstruct(x)));
p.addParameter('timestamp',char(datetime('now')),@ischar);
p.addParameter('hostname',char(java.lang.System.getProperty('user.name')),@ischar);
p.addParameter('username',char(java.net.InetAddress.getLocalHost.getHostName),@ischar);

%parse
p.parse(grayVideoName, glintFileName, varargin{:})


%% read video file
% load pupilPerimeter
videoInObj = VideoReader(grayVideoName);

% get number of frames
nFrames = floor(videoInObj.Duration*videoInObj.FrameRate);

% %% Initialize glint struct
% 
% % main glint params
% glintData.X = nan(nFrames,1);
% glintData.Y = nan(nFrames,1);
% glintData.size = nan(nFrames,1);
% 
% % glint fit params
% glintData.implicitEllipseParams = nan(nFrames,6);
% glintData.explicitEllipseParams= nan(nFrames,5);
% glintData.distanceErrorMetric= nan(nFrames,1);
% 
% % glint mask params
% glintData.circleRad = nan(nFrames,1);
% glintData.circleX = nan(nFrames,1);
% glintData.circleY = nan(nFrames,1);
% glintData.circleStrength = nan(nFrames,1);
% 
% % glint flags
% glintData.ellipseFittingError = nan(nFrames,1);


%% Track the glint

% alert the user
if strcmp(p.Results.verbosity,'full')
    tic
    fprintf(['Tracking the glint. Started ' char(datetime('now')) '\n']);
    fprintf('| 0                      50                   100%% |\n');
    fprintf('.');
end

%loop through frames
for ii = 1:nFrames
    % increment the progress bar
    if strcmp(p.Results.verbosity,'full') && mod(ii,round(nFrames/50))==0
        fprintf('.');
    end
    
    % Get the frame
    thisFrame = readFrame(videoInObj);
    
    % adjust gamma for this frame
    thisFrame = imadjust(thisFrame,[],[],p.Results.gammaCorrection);
    
    thisFrame = rgb2gray (thisFrame);
    % track with circles (using the default options)
    [~,~,~, gCenters, gRadii,gMetric, pupilRange, glintRange] = circleFit(thisFrame, ...
        p.Results.pupilCircleThresh, ...
        p.Results.glintCircleThresh, ...
        p.Results.pupilRange, ...
        p.Results.glintRange);
    
    % if a glint was present, refine the location
    if ~isempty(gCenters)
        % getGlintPerimeter
        [binG] = getGlintPerimeter (thisFrame, gCenters, gRadii, p.Results.glintEllipseThresh);
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
            glintData.X(ii)= gCenters(1,1);
            glintData.Y(ii) = gCenters(1,2);
            glintData.size(ii) = gRadii(1);
            glintData.circleStrength(ii) = gMetric(1);
            glintData.ellipseFittingError(ii) = 1;
            clear ME
        end
        
        % store results
        if exist ('Eg','var')
            if ~isempty (Eg) && isreal(Egi)
                glintData.X(ii) = Eg(2);
                glintData.Y(ii) = Eg(1);
                glintData.circleStrength(ii) = gMetric(1);
                glintData.implicitEllipseParams(ii,:) = Egi';
                glintData.explicitEllipseParams(ii,:) = Eg';
                glintData.distanceErrorMetric(ii) = gdistanceErrorMetric;
                % circle params for glint
                glintData.circleStrength(ii) = gMetric(1);
                glintData.circleRad(ii) = gRadii(1);
                glintData.circleX(ii) = gCenters(1,1);
                glintData.circleY(ii) = gCenters(1,2);
            end
            clear Eg Egi errors
        else
            glintData.X(ii)= gCenters(1,1);
            glintData.Y(ii) = gCenters(1,2);
            glintData.size(ii) = gRadii(1);
            glintData.circleStrength(ii) = gMetric(1);
        end
    end
end

close all
clear inObj

% add a meta field with analysis details
glintData.meta = p.Results;

%% save out a mat file with the glint tracking data
save (glintFileName, 'glintData')

% report completion of analysis
if strcmp(p.Results.verbosity,'full')
    fprintf('\n');
    toc
    fprintf('\n');
end

end % function


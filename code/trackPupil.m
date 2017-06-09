function [pupil, glint, params] = trackPupil(params)

% Tracks the pupil using an input video file, write out an .avi video
% By default, the video is resized and cropped to have the same aspect
% ratio and size of the original LiveTrack video input.
%
% The default fitting is done with an ellipse, after applying a circular
% mask to the original video. A simple circular fit can be used for
% tracking, but the resulting data appears noisier and generally less
% accurate.
%
%   Usage:
%       [pupil,glint,params]       = trackPupil(params)
%
%   Required input:
%       params.inVideo      = '/path/to/inputFile';
%
%   Optional inputs:
%       params.outVideo     = '/path/to/outVideo.avi';
%       params.outMat       = '/path/to/outData.mat'
%
%   Outputs:
%       pupil.X             = X coordinate of pupil center (pixels)
%       pupil.Y             = Y coordinate of pupil center (pixels)
%       pupil.size          = radius of pupil (pixels)
%       pupil.circleStrength= strength (accuracy) of the circular mask
%       pupil.ellipse       = parameters of the tracked ellipse
%       glint.X             = X coordinate of glint center (pixels)
%       glint.Y             = Y coordinate of glint center (pixels)
%       glint.size          = radius of glint (pixels)
%       glint.circleStrength= strength (accuracy) of the circular mask
%       glint.ellipse       = parameters of the tracked ellipse
%       params              = struct with all params used for tracking
%
%   Defaults:
%       params.pupilFit     = 'ellipse';    % tracking strategy
%       params.pupilOnly    = 0;            % if 1, no glint is required
%       params.imageSize    = [486 720]/2;; % used to resize input video (for no resizing, just input the native image size [ Y X ] )
%       params.imageCrop    = [1 1 319 239] % used to crop the image
%
%       params.rangeAdjust  = 0.05;         % radius change (+/-) allowed from the previous frame for circular mask
%       params.circleThresh   = [0.05 0.999]; % grayscale threshold values for pupil and glint, respectively
%       params.pupilRange   = [10 70];      % initial pupil size range
%       params.glintRange   = [10 30];      % constant glint size range
%       params.glintOut     = 0.1;          % proportion outside of pupil glint is allowed to be. Higher = more outside
%       params.sensitivity  = 0.99;         % [0 1] - sensitivity for 'imfindcircles'. Higher = more circles found
%       params.dilateGlint  = 6;            % used to dialate glint. Higher = more dilation.
%
%       params.ellipseThresh= [0.95 0.9];   % used to threshold the masked image for ellipse tracking
%       params.maskBox      = [4 8];        % used to expand the circular mask for ellipse tracking
%
%   V 2.0 -- Written by Giulia Frazzetta Feb 2017 : general code reorganization,
%   changed default video size, added ellipse option, output eyetracking
%   params.
%   V 1.0 -- Written by Andrew S Bock Sep 2016

%% set defaults

% params to choose the tracking
if ~isfield(params,'pupilFit')
    params.pupilFit = 'bestPupilCut';
end
if ~isfield(params,'pupilOnly')
    params.pupilOnly = 0;
end

% params for image resizing and cropping
if ~isfield (params, 'keepOriginalSize')
    params.keepOriginalSize = 0;
end
if ~isfield(params,'imageSize')
    params.imageSize = [486 720]/2;
end
if ~isfield(params,'imageCrop')
    params.imageCrop = [1 1 319 239];
end

% params for circleFit (always needed)
if ~isfield(params,'rangeAdjust')
    params.rangeAdjust = 0.05;
end
if ~isfield(params,'circleThresh')
    params.circleThresh = [0.06 0.999];
end
if ~isfield(params,'pupilRange')
    params.pupilRange   = [20 90];
end
if ~isfield(params,'glintRange')
    params.glintRange   = [10 30];
end
if ~isfield(params,'glintOut')
    params.glintOut     = 0.1;
end
if ~isfield(params,'sensitivity')
    params.sensitivity  = 0.99;
end
if ~isfield(params,'dilateGlint')
    params.dilateGlint  = 5;
end


% params for ellipse fit
if ~isfield(params,'ellipseThresh')
    params.ellipseThresh   = [0.97 0.9];
end
if ~isfield(params,'maskBox')
    params.maskBox   = [4 30];
end
if ~isfield(params,'gammaCorrection')
    params.gammaCorrection   = 1;
end
if ~isfield(params,'overGlintCut')
    params.overGlintCut = 5;
end
if ~isfield(params,'cutEveryFrame')
    params.cutEveryFrame = 0;
end
if ~isfield(params,'fittingErrorThresh')
    params.fittingErrorThresh = 10;
end

%% Load video
disp('Loading video file...');
inObj                   = VideoReader(params.inVideo);
numFrames               = floor(inObj.Duration*inObj.FrameRate);
% option to overwrite numFrames (for quick testing)
if isfield(params,'forceNumFrames')
    numFrames = params.forceNumFrames;
end

% initialize gray image array
grayI                   = zeros([240 320 numFrames],'uint8');

disp('Converting video to standard format, may take a while...');
% Convert to gray, resize, crop to livetrack size
for i = 1:numFrames
    thisFrame           = readFrame(inObj);
    tmp                 = rgb2gray(thisFrame);
    if params.keepOriginalSize == 0
        tmp2 = imresize(tmp,params.imageSize);
        tmp = imcrop(tmp2,params.imageCrop);
    end
    grayI(:,:,i) = tmp;
end


if isfield(params,'outVideo')
    outObj              = VideoWriter(params.outVideo);
    outObj.FrameRate    = inObj.FrameRate;
    open(outObj);
end

clear RGB inObj

%% Initialize pupil and glint structures
% display main tracking parameters
disp('Starting tracking with the following parameters:');
disp('Track pupil only: ')
disp(params.pupilOnly)
disp('Fit method: ')
disp(params.pupilFit)
disp('Circle threshold: ')
disp(params.circleThresh)
disp('Ellipse threshold: ')
disp(params.ellipseThresh)

% initialize
switch params.pupilFit
    case 'circle'
        % Useful:
        %   figure;imshow(I);
        %   d = imdistline;
        %   % check size of pupil or glint
        %   delete(d);
        
        pupilRange = params.pupilRange;
        glintRange = params.glintRange;
        pupil.X = nan(numFrames,1);
        pupil.Y = nan(numFrames,1);
        pupil.size = nan(numFrames,1);
        pupil.circleStrength = nan(numFrames,1);
        glint.X = nan(numFrames,1);
        glint.Y = nan(numFrames,1);
        glint.size = nan(numFrames,1);
        glint.circleStrength = nan(numFrames,1);
        
    case 'fixedPupilCut'
        % main pupil params
        pupil.X = nan(numFrames,1);
        pupil.Y = nan(numFrames,1);
        pupil.size = nan(numFrames,1);
        
        % full fit params
        pupil.implicitEllipseParams = nan(numFrames,6);
        pupil.explicitEllipseParams= nan(numFrames,5);
        pupil.distanceErrorMetric= nan(numFrames,1);
        pupil.pixelsOnPerimeter = nan(numFrames,1);
        
        % cut fit params
        pupil.cutPixels = nan(numFrames,1);
        pupil.cutImplicitEllipseParams = nan(numFrames,6);
        pupil.cutExplicitEllipseParams= nan(numFrames,5);
        pupil.cutDistanceErrorMetric= nan(numFrames,1);
        
        % pupil mask params
        pupilRange = params.pupilRange;
        pupil.circleRad = nan(numFrames,1);
        pupil.circleX = nan(numFrames,1);
        pupil.circleY = nan(numFrames,1);
        pupil.circleStrength = nan(numFrames,1);
        % structuring element for pupil mask size
        sep = strel('rectangle',params.maskBox);
        
        % pupil flags
        pupil.flags.fittingFailure = nan(numFrames,1);
        pupil.flags.noGlint = nan(numFrames,1);
        pupil.flags.noPupil = nan(numFrames,1);
        pupil.flags.cutPupil = nan(numFrames,1);
        pupil.flags.highCutError = nan(numFrames,1);
        
        % main glint params
        glint.X = nan(numFrames,1);
        glint.Y = nan(numFrames,1);
        glint.size = nan(numFrames,1);
        
        % glint fit params
        glint.implicitEllipseParams = nan(numFrames,6);
        glint.explicitEllipseParams= nan(numFrames,5);
        glint.distanceErrorMetric= nan(numFrames,1);
        
        % glint mask params
        glintRange = params.glintRange;
        glint.circleRad = nan(numFrames,1);
        glint.circleX = nan(numFrames,1);
        glint.circleY = nan(numFrames,1);
        glint.circleStrength = nan(numFrames,1);
        
        % glint flags
        glint.flags.fittingFailure = nan(numFrames,1);
        
    case 'bestPupilCut'
        
        % define cuts
        cuts = {...
            'full' ...
            'cut75up' ... % cuts 75% of the pixels over the glint
            'cut75right'... % cuts 75% of the pixels right of the glint 
            'cutDiagonal4075' ... 
            };
        
        % best pupil params
        pupil.X = nan(numFrames,1);
        pupil.Y = nan(numFrames,1);
        pupil.size = nan(numFrames,1);
        pupil.bestCut = cell(numFrames,1);
        pupil.bestError = nan(numFrames,1);
        pupil.fullPerimeterLength = nan(numFrames,1);
        pupil.area = nan(numFrames,1);
        
        % all fitting params
        for cc = 1 : length(cuts)
            % perimeter params
            pupil.(cuts{cc}).cutPixels = nan(numFrames,1);
            pupil.(cuts{cc}).flags.cutError = nan(numFrames,1);
            % ellipse params
            pupil.(cuts{cc}).implicitEllipseParams = nan(numFrames,6);
            pupil.(cuts{cc}).explicitEllipseParams= nan(numFrames,5);
            % area params
            pupil.(cuts{cc}).area = nan(numFrames,1);
            pupil.(cuts{cc}).varAreaPerSec = nan(numFrames,1);                   
            % error params
            pupil.(cuts{cc}).distanceError= nan(numFrames,1);
            pupil.(cuts{cc}).axRatio = nan(numFrames,1);
            pupil.(cuts{cc}).areaVariationError = nan(numFrames,1);
            pupil.(cuts{cc}).circularityError = nan(numFrames,1);
            % eccentricity
%             pupil.(cuts{cc}).Yeccentricity = nan(numFrames,1);
%             pupil.(cuts{cc}).Xeccentricity = nan(numFrames,1);
        end
        
        % pupil flags
        pupil.flags.fittingFailure = nan(numFrames,1);
        pupil.flags.noGlint = nan(numFrames,1);
        pupil.flags.noPupil = nan(numFrames,1);
        
        % pupil mask params
        pupilRange = params.pupilRange;
        pupil.circleRad = nan(numFrames,1);
        pupil.circleX = nan(numFrames,1);
        pupil.circleY = nan(numFrames,1);
        pupil.circleStrength = nan(numFrames,1);
        % structuring element for pupil mask size
        sep = strel('rectangle',params.maskBox);
        
        % main glint params
        glint.X = nan(numFrames,1);
        glint.Y = nan(numFrames,1);
        glint.size = nan(numFrames,1);
        
        % glint fit params
        glint.implicitEllipseParams = nan(numFrames,6);
        glint.explicitEllipseParams= nan(numFrames,5);
        glint.distanceErrorMetric= nan(numFrames,1);
        
        % glint mask params
        glintRange = params.glintRange;
        glint.circleRad = nan(numFrames,1);
        glint.circleX = nan(numFrames,1);
        glint.circleY = nan(numFrames,1);
        glint.circleStrength = nan(numFrames,1);
        
        % glint flags
        glint.flags.fittingFailure = nan(numFrames,1);
        
        
        
end

%% Track
progBar = ProgressBar(numFrames,'tracking pupil...');
if isfield(params,'outVideo')
    ih = figure;
end

switch params.pupilFit
    
    case 'circle'
        for i = 1:numFrames
            % Get the frame
            I = squeeze(grayI(:,:,i));
            % Show the frame
            if isfield(params,'outVideo')
                imshow(I);
            end
            
            % track with circles
            [pCenters, pRadii,pMetric, gCenters, gRadii,gMetric, pupilRange, glintRange] = circleFit(I,params,pupilRange,glintRange);
            
            % visualize tracking results
            switch params.pupilOnly
                case 0
                    if ~isempty(pCenters) && ~isempty(gCenters)
                        pupil.X(i) = pCenters(1,1);
                        pupil.Y(i) = pCenters(1,2);
                        pupil.size(i) = pRadii(1);
                        pupil.circleStrength(i) = pMetric(1);
                        glint.X(i)= gCenters(1,1);
                        glint.Y(i) = gCenters(1,2);
                        glint.size(i) = gRadii(1);
                        glint.circleStrength(i) = gMetric(1);
                        if isfield(params,'outVideo')
                            viscircles(pCenters(1,:),pRadii(1),'Color','r');
                            hold on
                            plot(gCenters(1,1),gCenters(1,2),'+b');
                            hold off
                        end
                    end
                case 1
                    if ~isempty(pCenters)
                        pupil.X(i)      = pCenters(1,1);
                        pupil.Y(i)      = pCenters(1,2);
                        pupil.size(i)   = pRadii(1);
                        pupil.strength(i)  = pMetric(1);
                        if isfield(params,'outVideo')
                            viscircles(pCenters(1,:),pRadii(1),'Color','r');
                        end
                    end
            end
            
            % save single frame
            if isfield(params,'outVideo')
                frame                   = getframe(ih);
                writeVideo(outObj,frame);
            end
            if ~mod(i,10);progBar(i);end;
        end
        
    case 'fixedPupilCut'
        for i = 1:numFrames
            % Get the frame
            I = squeeze(grayI(:,:,i));
            % adjust gamma for this frame
            I = imadjust(I,[],[],params.gammaCorrection);
%             I = imadjust(I,[.22 .4],[]);
            % Show the frame
            if isfield(params,'outVideo')
                imshow(I);
            end
            
            % track with circles
            [pCenters, pRadii,pMetric, gCenters, gRadii,gMetric, pupilRange, glintRange] = circleFit(I,params,pupilRange,glintRange);
            
            switch params.pupilOnly
                
                % track pupil only
                case 1
                    if isempty(pCenters)
                        % save frame
                        if isfield(params,'outVideo')
                            frame   = getframe(ih);
                            writeVideo(outObj,frame);
                        end
                        if ~mod(i,10);progBar(i);end;
                        continue
                    else
                        % get pupil perimeter
                        [binP] = getPupilPerimeter(I,pCenters,pRadii, sep, params);
                        try
                            % Fit ellipse to pupil
                            [Xp, Yp] = ind2sub(size(binP),find(binP));
                            Epi = ellipsefit_direct(Xp,Yp);
                            Ep = ellipse_im2ex(Epi);
                            
                            % get errorMetric
                            [~,d,~,~] = ellipse_distance(Xp, Yp, Epi);
                            distanceErrorMetric = nanmedian(sqrt(sum(d.^2)));
                        catch ME
                        end
                        if  exist ('ME', 'var')
                            % save frame
                            if isfield(params,'outVideo')
                                frame   = getframe(ih);
                                writeVideo(outObj,frame);
                            end
                            if ~mod(i,10);progBar(i);end;
                            clear ME
                            continue
                        end
                            
                        % store results
                        if ~isempty(Ep)
                            pupil.X(i) = Ep(2);
                            pupil.Y(i) = Ep(1);
                            pupil.size(i) = Ep(3); % " bigger radius"
                            % ellipse params
                            pupil.implicitEllipseParams(i,:) = Epi';
                            pupil.explicitEllipseParams(i,:) = Ep';
                            pupil.distanceErrorMetric(i) = distanceErrorMetric;
                            % circle patch params
                            pupil.circleStrength(i) = pMetric(1);
                            pupil.circleRad(i) = pRadii(1);
                            pupil.circleX(i) = pCenters(1,1);
                            pupil.circleY(i) = pCenters(1,2);
                            
                            % plot results
                             if ~isempty(Epi) && Ep(1) > 0
                                a = num2str(Epi(1));
                                b = num2str(Epi(2));
                                c = num2str(Epi(3));
                                d = num2str(Epi(4));
                                e = num2str(Epi(5));
                                f = num2str(Epi(6));
                                
                                % note that X and Y indices need to be swapped!
                                eqt= ['(',a, ')*y^2 + (',b,')*x*y + (',c,')*x^2 + (',d,')*y+ (',e,')*x + (',f,')'];
                                
                                if isfield(params,'outVideo')
                                    hold on
                                    h= ezplot(eqt,[1, 240, 1, 320]);
                                    % set color according to type of tracking
                                    set (h, 'Color', 'green')
                                    hold off
                                end
                            end % plot results
                             % save frame
                            if isfield(params,'outVideo')
                                frame   = getframe(ih);
                                writeVideo(outObj,frame);
                            end
                            if ~mod(i,10);progBar(i);end;
                        else
                            % circle params
                            pupil.circleStrength(i) = pMetric(1);
                            pupil.circleRad(i) = pRadii(1);
                            pupil.circleX(i) = pCenters(1,1);
                            pupil.circleY(i) = pCenters(1,2);
                            
                            % save frame
                            if isfield(params,'outVideo')
                                frame   = getframe(ih);
                                writeVideo(outObj,frame);
                            end
                            if ~mod(i,10);progBar(i);end;
                            continue
                        end
                    end
                    
                    % track pupil and glint
                case 0
                    if isempty(gCenters) && isempty(pCenters)
                        pupil.flags.noPupil(i) = 1;
                        pupil.flags.noGlint(i) = 1;
                        % save frame
                        if isfield(params,'outVideo')
                            frame   = getframe(ih);
                            writeVideo(outObj,frame);
                        end
                        if ~mod(i,10);progBar(i);end;
                        continue
                    elseif isempty(pCenters) && ~isempty(gCenters)
                        pupil.flags.noPupil(i) = 1;
                        pupil.flags.noGlint(i) = 0;
                        % circle params for glint
                        glint.circleStrength(i) = gMetric(1);
                        glint.circleRad(i) = gRadii(1);
                        glint.circleX(i) = gCenters(1,1);
                        glint.circleY(i) = gCenters(1,2);
                        % save frame
                        if isfield(params,'outVideo')
                            frame   = getframe(ih);
                            writeVideo(outObj,frame);
                        end
                        if ~mod(i,10);progBar(i);end;
                        continue
                        
                    elseif ~isempty(pCenters) && isempty(gCenters)
                        pupil.flags.noPupil(i) = 0;
                        pupil.flags.noGlint(i) = 1;
                        % circle params for pupil
                        pupil.circleStrength(i) = pMetric(1);
                        pupil.circleRad(i) = pRadii(1);
                        pupil.circleX(i) = pCenters(1,1);
                        pupil.circleY(i) = pCenters(1,2);
                        % save frame
                        if isfield(params,'outVideo')
                            frame   = getframe(ih);
                            writeVideo(outObj,frame);
                        end
                        if ~mod(i,10);progBar(i);end;
                        continue
                        
                    elseif ~isempty(pCenters) && ~isempty(gCenters)
                        pupil.flags.noPupil(i) = 0;
                        pupil.flags.noGlint(i) = 0;
                        
                        % Track the glint
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
                            glint.flags.fittingError(i) = 1;
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
                        
                        
                        % Track the pupil
                        % get pupil perimeter
                        [binP] = getPupilPerimeter(I,pCenters,pRadii, sep, params);
                        [Xp, Yp] = ind2sub(size(binP),find(binP));
                        params.pixelsOnPerimeter(i) = length(Xp);
                        
                        if params.cutEveryFrame % cut the perimeter for every frame, regardless of error.
                            pupil.flags.cutPupil(i) = 1;
                            underGlint = find (Xp > gCenters(1,2) - params.overGlintCut );
                            params.cutPixels(i) = length(underGlint);
                            if ~isempty(underGlint)
                                binPcut = zeros(size(I));
                                binPcut(sub2ind(size(binP),Xp(underGlint),Yp(underGlint))) = 1;
                                %                               imshow(binPcut)
                                binP = binPcut;
                                [Xp, Yp] = ind2sub(size(binP),find(binP));
                            end
                        end
                        
                        % Fit ellipse to pupil
                        try
                            Epi = ellipsefit_direct(Xp,Yp);
                            Ep = ellipse_im2ex(Epi);
                            
                            % get errorMetric
                            [~,d,~,~] = ellipse_distance(Xp, Yp, Epi);
                            distanceErrorMetric = nanmedian(sqrt(sum(d.^2)));
                        catch ME
                        end
                        if  exist ('ME', 'var')
                            % circle params
                            pupil.circleStrength(i) = pMetric(1);
                            pupil.circleRad(i) = pRadii(1);
                            pupil.circleX(i) = pCenters(1,1);
                            pupil.circleY(i) = pCenters(1,2);
                            pupil.flags.fittingError(i) = 1;
                            % save frame
                            if isfield(params,'outVideo')
                                frame   = getframe(ih);
                                writeVideo(outObj,frame);
                            end
                            if ~mod(i,10);progBar(i);end;
                            clear ME
                            continue
                        end
                        
                        % check the error of the fitting
                        if params.cutEveryFrame
                            % store results
                            if exist ('Ep','var')
                                if ~isempty(Ep) && isreal(Epi)
                                    pupil.X(i) = Ep(2);
                                    pupil.Y(i) = Ep(1);
                                    pupil.size(i) = Ep(3); % "radius"
                                    % ellipse params
                                    pupil.cutImplicitEllipseParams(i,:) = Epi';
                                    pupil.cutExplicitEllipseParams(i,:) = Ep';
                                    pupil.cutDistanceErrorMetric(i) = distanceErrorMetric;
                                    % circle params
                                    pupil.circleStrength(i) = pMetric(1);
                                    pupil.circleRad(i) = pRadii(1);
                                    pupil.circleX(i) = pCenters(1,1);
                                    pupil.circleY(i) = pCenters(1,2);
                                end
                                if distanceErrorMetric > params.fittingErrorThresh
                                    pupil.flags.highCutError(i) = 1;
                                else
                                    pupil.flags.highCutError(i) = 0;
                                end
                            else
                                % circle params
                                pupil.circleStrength(i) = pMetric(1);
                                pupil.circleRad(i) = pRadii(1);
                                pupil.circleX(i) = pCenters(1,1);
                                pupil.circleY(i) = pCenters(1,2);
                                % save frame
                                if isfield(params,'outVideo')
                                    frame   = getframe(ih);
                                    writeVideo(outObj,frame);
                                end
                                if ~mod(i,10);progBar(i);end;
                                continue
                            end
                        elseif ~params.cutEveryFrame
                            % store results
                            if exist ('Ep','var')
                                if ~isempty(Ep) && isreal(Epi)
                                    pupil.X(i) = Ep(2);
                                    pupil.Y(i) = Ep(1);
                                    pupil.size(i) = Ep(3); % "radius"
                                    % ellipse params
                                    pupil.implicitEllipseParams(i,:) = Epi';
                                    pupil.explicitEllipseParams(i,:) = Ep';
                                    pupil.distanceErrorMetric(i) = distanceErrorMetric;
                                    % circle params
                                    pupil.circleStrength(i) = pMetric(1);
                                    pupil.circleRad(i) = pRadii(1);
                                    pupil.circleX(i) = pCenters(1,1);
                                    pupil.circleY(i) = pCenters(1,2);
                                end
                                pupil.flags.cutPupil(i) = 0;
                                pupil.flags.highCutError(i) = 0;
                            else
                                % circle params
                                pupil.circleStrength(i) = pMetric(1);
                                pupil.circleRad(i) = pRadii(1);
                                pupil.circleX(i) = pCenters(1,1);
                                pupil.circleY(i) = pCenters(1,2);
                                % save frame
                                if isfield(params,'outVideo')
                                    frame   = getframe(ih);
                                    writeVideo(outObj,frame);
                                end
                                if ~mod(i,10);progBar(i);end;
                                continue
                            end
                            
                            if distanceErrorMetric > params.fittingErrorThresh
                                % clear previous ellipse variables
                                clear Ep Epi
                                % cut this frame
                                pupil.flags.cutPupil(i) = 1;
                                underGlint = find (Xp > gCenters(1,2) - params.overGlintCut );
                                params.cutPixels(i) = length(underGlint);
                                if ~isempty(underGlint)
                                    binPcut = zeros(size(I));
                                    binPcut(sub2ind(size(binP),Xp(underGlint),Yp(underGlint))) = 1;
                                    %                               imshow(binPcut)
                                    binP = binPcut;
                                    [Xp, Yp] = ind2sub(size(binP),find(binP));
                                end
                                
                                % fit again
                                try
                                    Epi = ellipsefit_direct(Xp,Yp);
                                    Ep = ellipse_im2ex(Epi);
                                    
                                    % get errorMetric
                                    [~,d,~,~] = ellipse_distance(Xp, Yp, Epi);
                                    cutdistanceErrorMetric = nanmedian(sqrt(sum(d.^2)));
                                catch ME
                                end
                                if  exist ('ME', 'var')
                                    pupil.flags.fittingError(i) = 1;
                                    % save frame
                                    if isfield(params,'outVideo')
                                        frame   = getframe(ih);
                                        writeVideo(outObj,frame);
                                    end
                                    if ~mod(i,10);progBar(i);end;
                                    clear ME
                                    continue
                                end
                                % store results
                                if exist ('Ep','var')
                                    if ~isempty(Ep) && isreal(Epi)
                                        pupil.X(i) = Ep(2);
                                        pupil.Y(i) = Ep(1);
                                        pupil.size(i) = Ep(3); % "radius"
                                        % ellipse params
                                        pupil.cutImplicitEllipseParams(i,:) = Epi';
                                        pupil.cutExplicitEllipseParams(i,:) = Ep';
                                        pupil.cutDistanceErrorMetric(i) = cutdistanceErrorMetric;
                                    end
                                    if cutdistanceErrorMetric > params.fittingErrorThresh
                                        pupil.flags.highCutError(i) = 1;
                                    else
                                        pupil.flags.highCutError(i) = 0;
                                    end
                                end
                            end
                            
                            % plot results
                            if ~isempty(Epi) && Ep(1) > 0
                                a = num2str(Epi(1));
                                b = num2str(Epi(2));
                                c = num2str(Epi(3));
                                d = num2str(Epi(4));
                                e = num2str(Epi(5));
                                f = num2str(Epi(6));
                                
                                % note that X and Y indices need to be swapped!
                                eqt= ['(',a, ')*y^2 + (',b,')*x*y + (',c,')*x^2 + (',d,')*y+ (',e,')*x + (',f,')'];
                                
                                if isfield(params,'outVideo')
                                    hold on
                                    h= ezplot(eqt,[1, 240, 1, 320]);
                                    % set color according to type of tracking
                                    if pupil.flags.cutPupil(i) == 0
                                        set (h, 'Color', 'green')
                                    elseif pupil.flags.cutPupil(i) == 1 && pupil.flags.highCutError(i) == 0
                                        set (h, 'Color', 'yellow')
                                    elseif pupil.flags.cutPupil(i) == 1 && pupil.flags.highCutError(i) == 1
                                        set (h, 'Color', 'red')
                                    end
                                    if ~params.pupilOnly && ~isnan(glint.X(i))
                                        hold on
                                        plot(glint.X(i),glint.Y(i),'+b');
                                    end
                                    hold off
                                end
                            end % plot results
                        end
                    end
                    
                    % save frame
                    if isfield(params,'outVideo')
                        frame   = getframe(ih);
                        writeVideo(outObj,frame);
                    end
                    if ~mod(i,10);progBar(i);end;
                    clear Eg Egi Ep Epi
            end
        end
        % save full video
        if isfield(params,'outVideo')
            close(ih);
            close(outObj);
        end
        
        % save tracked values to output matrix
        if isfield(params,'outMat')
            save(params.outMat,'pupil','glint');
        end
        
        
        
        % -------- NEW CASE
        
    case 'bestPupilCut'
        for i = 1:numFrames
            % Get the frame
            I = squeeze(grayI(:,:,i));
            
            % adjust gamma for this frame
            I = imadjust(I,[],[],params.gammaCorrection);
           
            % Show the frame
            if isfield(params,'outVideo')
                imshow(I);
            end
            
            % track with circles
            [pCenters, pRadii,pMetric, gCenters, gRadii,gMetric, pupilRange, glintRange] = circleFit(I,params,pupilRange,glintRange);
            
            switch params.pupilOnly
                
                % track pupil only
                case 1
                    if isempty(pCenters)
                        % save frame
                        if isfield(params,'outVideo')
                            frame   = getframe(ih);
                            writeVideo(outObj,frame);
                        end
                        if ~mod(i,10);progBar(i);end;
                        continue
                    else
                        % get pupil perimeter
                        [binP] = getPupilPerimeter(I,pCenters,pRadii, sep, params);
                        
                        % Fit ellipse to pupil
                        [Xp, Yp] = ind2sub(size(binP),find(binP));
                        Epi = ellipsefit_direct(Xp,Yp);
                        Ep = ellipse_im2ex(Epi);
                        
                        % get errorMetric
                        [~,d,~,~] = ellipse_distance(Xp, Yp, Epi);
                        distanceErrorMetric = nanmedian(sqrt(sum(d.^2)));
                        
                        % store results
                        if ~isempty(Ep)
                            pupil.X(i) = Ep(2);
                            pupil.Y(i) = Ep(1);
                            pupil.size(i) = Ep(3); % " bigger radius"
                            % ellipse params
                            pupil.implicitEllipseParams(i,:) = Epi';
                            pupil.explicitEllipseParams(i,:) = Ep';
                            pupil.distanceErrorMetric(i) = distanceErrorMetric;
                            % circle patch params
                            pupil.circleStrength(i) = pMetric(1);
                            pupil.circleRad(i) = pRadii(1);
                            pupil.circleX(i) = pCenters(1,1);
                            pupil.circleY(i) = pCenters(1,2);
                        else
                            % circle params
                            pupil.circleStrength(i) = pMetric(1);
                            pupil.circleRad(i) = pRadii(1);
                            pupil.circleX(i) = pCenters(1,1);
                            pupil.circleY(i) = pCenters(1,2);
                            
                            % save frame
                            if isfield(params,'outVideo')
                                frame   = getframe(ih);
                                writeVideo(outObj,frame);
                            end
                            if ~mod(i,10);progBar(i);end;
                            continue
                        end
                    end
                    
                    % track pupil and glint
                case 0
                    if isempty(gCenters) && isempty(pCenters)
                        pupil.flags.noPupil(i) = 1;
                        pupil.flags.noGlint(i) = 1;
                        % save frame
                        if isfield(params,'outVideo')
                            frame   = getframe(ih);
                            writeVideo(outObj,frame);
                        end
                        if ~mod(i,10);progBar(i);end;
                        continue
                    elseif isempty(pCenters) && ~isempty(gCenters)
                        pupil.flags.noPupil(i) = 1;
                        pupil.flags.noGlint(i) = 0;
                        % circle params for glint
                        glint.circleStrength(i) = gMetric(1);
                        glint.circleRad(i) = gRadii(1);
                        glint.circleX(i) = gCenters(1,1);
                        glint.circleY(i) = gCenters(1,2);
                        % save frame
                        if isfield(params,'outVideo')
                            frame   = getframe(ih);
                            writeVideo(outObj,frame);
                        end
                        if ~mod(i,10);progBar(i);end;
                        continue
                        
                    elseif ~isempty(pCenters) && isempty(gCenters)
                        pupil.flags.noPupil(i) = 0;
                        pupil.flags.noGlint(i) = 1;
                        % circle params for pupil
                        pupil.circleStrength(i) = pMetric(1);
                        pupil.circleRad(i) = pRadii(1);
                        pupil.circleX(i) = pCenters(1,1);
                        pupil.circleY(i) = pCenters(1,2);
                        % save frame
                        if isfield(params,'outVideo')
                            frame   = getframe(ih);
                            writeVideo(outObj,frame);
                        end
                        if ~mod(i,10);progBar(i);end;
                        continue
                        
                    elseif ~isempty(pCenters) && ~isempty(gCenters)
                        pupil.flags.noPupil(i) = 0;
                        pupil.flags.noGlint(i) = 0;
                        
                        % Track the glint
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
                            glint.flags.fittingError(i) = 1;
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
                        
                        
                        % ------ NEW STUFF -------
                        % Track the pupil
                        % get pupil perimeter
                        [binP] = getPupilPerimeter(I,pCenters,pRadii, sep, params);
                        [Xp, Yp] = ind2sub(size(binP),find(binP));
                        
                        % store full perimeter length
                        pupil.fullPerimeterLength(i) = length(Xp);
                        
                        % find pixels under and over glint
                        underGlint = find (Xp > gCenters(1,2));  %%% MUST CHECK IF EMPTY
                        overGlint = find (Xp < gCenters(1,2)); %%% MUST CHECK IF EMPTY
                        leftGlint = find (Yp < gCenters(1,1));
                        rightGlint = find (Yp > gCenters(1,1));
                        
                        
                        for cc = 1 : length(cuts) % loop through cuts
                            cut = cuts{cc};
                            switch cut
                                case 'full' % do not cut
                                    % perimeter params
                                    Xc = Xp;
                                    Yc = Yp;
                                    pupil.(cuts{cc}).cutPixels(i) = 0;
                                case 'cut75up'
                                    % cut out top 75% of the overglint
                                    cutout = round(length(overGlint)/100 * 25);
                                    pupil.(cuts{cc}).cutPixels(i) = cutout;
                                    if cutout > 0
                                        [~, idx] = sort(overGlint);
                                        % get the cut perimeter
                                        binPcut = zeros(size(I));
                                        binPcut(sub2ind(size(binP),Xp(underGlint),Yp(underGlint))) = 1;
                                        binPcut(sub2ind(size(binP),Xp(overGlint),Yp(overGlint))) = 1;
                                        binPcut(sub2ind(size(binP),Xp(overGlint(idx(1+round(cutout/2):end - round(cutout/2)))),Yp(overGlint(idx(1+round(cutout/2):end - round(cutout/2)))))) = 0;
                                        % remove small objects 
                                        CC = bwconncomp(binPcut);
                                        numPixels = cellfun(@numel,CC.PixelIdxList);
                                        [biggest,idx] = max(numPixels);
                                        if idx >0
                                            binPcut = zeros(size(I));
                                            binPcut(CC.PixelIdxList{idx}) = 1;
                                        end
                                        [Xc, Yc] = ind2sub(size(binPcut),find(binPcut));
                                        pupil.(cuts{cc}).flags.cutError(i) = 0;
                                    else
                                        pupil.(cuts{cc}).flags.cutError(i) = 1;
                                        continue
                                    end
                                case 'cut75right'
                                    cutout = round(length(rightGlint)/100 * 25);
                                    pupil.(cuts{cc}).cutPixels(i) = cutout;
                                    if cutout > 0
                                        [~, idx] = sort(rightGlint);
                                        % get the cut perimeter
                                        binPcut = zeros(size(I));
                                        binPcut(sub2ind(size(binP),Xp(leftGlint),Yp(leftGlint))) = 1;
                                        binPcut(sub2ind(size(binP),Xp(rightGlint),Yp(rightGlint))) = 1;
                                        binPcut(sub2ind(size(binP),Xp(rightGlint(idx(cutout:end))),Yp(rightGlint(idx(cutout:end))))) = 0;
                                        % remove small objects 
                                        CC = bwconncomp(binPcut);
                                        numPixels = cellfun(@numel,CC.PixelIdxList);
                                        [biggest,idx] = max(numPixels);
                                        if idx >0
                                            binPcut = zeros(size(I));
                                            binPcut(CC.PixelIdxList{idx}) = 1;
                                        end
                                        [Xc, Yc] = ind2sub(size(binPcut),find(binPcut));
                                        pupil.(cuts{cc}).flags.cutError(i) = 0;
                                    else
                                        pupil.(cuts{cc}).flags.cutError(i) = 1;
                                        continue
                                    end
                                    case 'cutDiagonal4075'
                                    cutoutU = round(length(overGlint)/100 * 60);
                                    cutoutR = round(length(rightGlint)/100 * 25);
                                    
                                    pupil.(cuts{cc}).cutPixels(i) = cutoutR + cutoutU;
                                    try
                                       binPcut = pupilDiagonalCut (I, cutoutU, cutoutR, binP, underGlint, overGlint, leftGlint, rightGlint);
                                        
                                        [Xc, Yc] = ind2sub(size(binPcut),find(binPcut));
                                        pupil.(cuts{cc}).flags.cutError(i) = 0;
                                    catch ME
                                    end
                                    if  exist ('ME', 'var')
                                        pupil.(cuts{cc}).flags.cutError(i) = 1;
                                        clear ME
                                        continue
                                    end
                            end % switch cut
                            
                            % do the fitting on the surviving perimeter pixels
                            try
                                Epi = ellipsefit_direct(Xc,Yc);
                                Ep = ellipse_im2ex(Epi);
                                
                                % store error parameters
                                [~,d,~,~] = ellipse_distance(Xc, Yc, Epi);
                                pupil.(cuts{cc}).distanceError(i) = nanmedian(sqrt(sum(d.^2)));
                                pupil.(cuts{cc}).area(i) = pi * Ep(3)* Ep(4);
                                pupil.(cuts{cc}).axRatio(i) = Ep(3)./Ep(4);

                                pupil.(cuts{cc}).circularityError(i) = 1+ 1 ./ (1+exp( -(pupil.(cuts{cc}).axRatio(i)*20-28) ))*9;
                                
                                previousFrames = 6;
                                if i > previousFrames
                                    % get last frames best errors
                                    bestErrors = pupil.bestError(i-previousFrames:i-1);
                                    prevAreas = pupil.area(i-previousFrames:i-1);
                                    bestIDX =  find(bestErrors < 20);
                                    if ~isempty (bestIDX)
                                        prevArea = nanmedian(prevAreas(bestIDX));
                                        pupil.(cuts{cc}).varAreaPerSec(i) = ((pupil.(cuts{cc}).area(i) - prevArea) /prevArea) *10 * 100;
                                        pupil.(cuts{cc}).areaVariationError(i) = ( 1 ./ (1+exp( -(abs(pupil.(cuts{cc}).varAreaPerSec(i))*.25-12.5) )))*7+1;
                                    elseif isempty (bestIDX) && ~isnan(pupil.(cuts{cc}).areaVariationError(i-1))
                                        pupil.(cuts{cc}).areaVariationError(i) = pupil.(cuts{cc}).areaVariationError(i-1);
                                    else
                                        pupil.(cuts{cc}).areaVariationError(i) = 1;
                                    end
                                    
                                else
                                   pupil.(cuts{cc}).areaVariationError(i) = 1;
                                end

                            catch ME
                            end
                            if  exist ('ME', 'var')
                                pupil.(cuts{cc}).flags.fittingError(i) = 1;
                                clear ME
                                continue
                            end
                            % store results
                            if exist ('Ep','var')
                                if ~isempty(Ep) && isreal(Epi)
                                    % ellipse params
                                    pupil.(cuts{cc}).implicitEllipseParams(i,:) = Epi';
                                    pupil.(cuts{cc}).explicitEllipseParams(i,:) = Ep';
                                else
                                    continue
                                end
                            else
                                pupil.(cuts{cc}).flags.fittingError(i) = 1;
                            end
                        end
                        
                        
                        % store other results
                        % circle params
                        pupil.circleStrength(i) = pMetric(1);
                        pupil.circleRad(i) = pRadii(1);
                        pupil.circleX(i) = pCenters(1,1);
                        pupil.circleY(i) = pCenters(1,2);
                        
                        % if the fit was good on the full perimeter just
                        % get that, otherwise compare to the cuts
                        if pupil.full.distanceError(i) < 12
                            eIDX = 1;
                            pupil.bestError(i) = pupil.(cuts{1}).distanceError(i) .* pupil.(cuts{1}).circularityError(i) .* pupil.(cuts{1}).areaVariationError(i);
                        else
                            % find best error
                            for ee = 1: length(cuts)
                                try
                                    errors(ee) = pupil.(cuts{ee}).distanceError(i) .* pupil.(cuts{ee}).circularityError(i) .* pupil.(cuts{ee}).areaVariationError(i);
                                catch ME
                                end
                                if exist ('ME', 'var')
                                    errors(ee) = NaN;
                                    clear ME
                                    continue
                                end
                            end
                            
                            if isempty(errors)
                                pupil.flags.fittingError(i) = 1;
                                % save frame
                                if isfield(params,'outVideo')
                                    frame   = getframe(ih);
                                    writeVideo(outObj,frame);
                                end
                                if ~mod(i,10);progBar(i);end;
                                continue
                            end
                            [pupil.bestError(i), eIDX] = nanmin (errors);
                        end
                        
                        pupil.bestCut(i) = {(cuts{eIDX})};
                        
                        % store results
                        pupil.X(i) = pupil.(cuts{eIDX}).explicitEllipseParams(i,2);
                        pupil.Y(i) = pupil.(cuts{eIDX}).explicitEllipseParams(i,1);
                        pupil.size(i) = pupil.(cuts{eIDX}).explicitEllipseParams(i,3);
                        pupil.area(i) = pupil.(cuts{eIDX}).area(i);
                        
                        % plot results
                        Epi = pupil.(cuts{eIDX}).implicitEllipseParams(i,:)' ;
                        Ep = pupil.(cuts{eIDX}).explicitEllipseParams(i,:)' ;
                        
                        if ~isempty(Epi) && Ep(1) > 0
                            a = num2str(Epi(1));
                            b = num2str(Epi(2));
                            c = num2str(Epi(3));
                            d = num2str(Epi(4));
                            e = num2str(Epi(5));
                            f = num2str(Epi(6));
                            
                            % note that X and Y indices need to be swapped!
                            eqt= ['(',a, ')*y^2 + (',b,')*x*y + (',c,')*x^2 + (',d,')*y+ (',e,')*x + (',f,')'];
                            
                            if isfield(params,'outVideo')
                                hold on
                                h= ezplot(eqt,[1, 240, 1, 320]);
                                % set color according to type of tracking
                                    if pupil.bestError(i) < 15
                                        set (h, 'Color', 'green')
                                    elseif pupil.bestError(i) < 50 && pupil.bestError(i) >= 15
                                        set (h, 'Color', 'yellow')
                                    elseif pupil.bestError(i) > 50
                                        set (h, 'Color', 'red')
                                    end
                                    if ~params.pupilOnly && ~isnan(glint.X(i))
                                        hold on
                                        plot(glint.X(i),glint.Y(i),'+b');
                                    end
                                    hold off
                            end
                        end
                        
                        % save frame
                        if isfield(params,'outVideo')
                            frame   = getframe(ih);
                            writeVideo(outObj,frame);
                        end
                        if ~mod(i,10);progBar(i);end;
                        clear Eg Egi Ep Epi errors
                    end
            end
        end
        
        % save video
        if isfield(params,'outVideo')
            close(ih);
            close(outObj);
        end
        
        % save tracked values to output matrix
        if isfield(params,'outMat')
            save(params.outMat,'pupil','glint');
        end
        
end






%% clean up
close all
clear
clc

%% define the following
[~, tmpName] = system('whoami');
userName = strtrim(tmpName);

dropboxDir = ['/Users/' userName '/Dropbox-Aguirre-Brainard-Lab'];
dropboxDir = '~/Desktop';

outDir = ['/Users/' userName '/Desktop'];
dropboxDir=outDir;

%% sample run
params.subjectName = 'TOME_3008';
params.sessionDate = '102116';
params.projectSubfolder = 'session1_restAndStructure';
params.runName = 'rfMRI_REST_AP_run01';

% params.subjectName = 'TOME_3013';
% params.sessionDate = '021717';
% params.projectSubfolder = 'session2_spatialstimuli';
% params.runName = 'tfMRI_FLASH_AP_run01';

%% set paths
videoPath = fullfile(dropboxDir,'TOME_processing',params.projectSubfolder,params.subjectName,params.sessionDate,'EyeTracking');
videoPath = dropboxDir;
params.inVideo = fullfile(videoPath,[params.runName '_60hz.avi']);
params.outVideo = fullfile(outDir,[params.runName '_perimeter.avi']);



%% tracking params
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
params.rangeAdjust = 0.05;
params.circleThresh = [0.085 0.999];
params.pupilRange   = [20 90];
params.glintRange   = [10 30];
params.glintOut     = 0.1;
params.sensitivity  = 0.99;
params.dilateGlint  = 5;
params.pupilOnly = 0;

% structuring element for pupil mask size
params.maskBox   = [4 30];
sep = strel('rectangle',params.maskBox);

% force number of frames
params.forceNumFrames = 1000;
params.ellipseThresh   = [0.9 0.9];
params.gammaCorrection = 1;


%% EXTRACTION OF THE PUPIL PERIMETER
% Load video
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
    ih = figure;
end
if isfield(params,'outVideo')
    outObj              = VideoWriter(params.outVideo);
    outObj.FrameRate    = inObj.FrameRate;
    open(outObj);
end

clear RGB inObj


%%%%%%%%%%%%%%%%%%%%%%%
% BEGIN GEOFF CODE HERE
%%%%%%%%%%%%%%%%%%%%%%%


% We will fit ellipses that are cast in "transparent" parameter form:
% center (cx,cy), its area (a), its eccentricity (e), and its angle of tilt
% (theta).

% Define the hard upper and lower boundaries.
% We can set an upper boundary on the eccentricity parameter. Eccentricity
% is related to ratio of the semimajor and semiminor axes, and can be
% calculated using:
%   eccentricity = axes2ecc(semimajor, semiminor)
% If we wish to prevent ellipses with an aspect ratio greater than 1.2 : 1,
% this gives us an eccentricity threshold of ~0.55.

lb = [0,  0,  1000,   0,  -0.5*pi];
ub = [240,320,10000,0.55, 0.5*pi];

% Define the initial prior values. Set to large SD values to cause the
% initial parameters to be dominated by the fits, and not the prior
pPriorMeanTransparent = [120,120,6000,0,0];
pPriorSDTransparent = [50,50,50,50,50];

% Define a prior window in units of samples
window=70;
windowSupport=1:1:window;
exponentialTauParam = round(window/4);

% Define a decaying exponential function that will be used to weight the
% contribution of prior values to the creation of the prior distribution
exponentialWeights=fliplr(exp(-1/exponentialTauParam*windowSupport));

% extract perimeter
for i = 1:numFrames
    % Get the frame
    I = squeeze(grayI(:,:,i));
    
    % adjust gamma for this frame
    I = imadjust(I,[],[],params.gammaCorrection);
    
    % track with circles
    pupilRange = params.pupilRange;
    glintRange = params.glintRange;
    [pCenters, pRadii,pMetric, gCenters, gRadii,gMetric, pupilRange, glintRange] = circleFit(I,params,pupilRange,glintRange);

    % get pupil perimeter
    [binP] = getPupilPerimeter(I,pCenters,pRadii, sep, params);

    % get the boundary points
    [Xc, Yc] = ind2sub(size(binP),find(binP));

    % fit an ellipse to the boundaary
    [pInitialFitTransparent, pFitSD, ~] = calcPupilLikelihood(Xc,Yc, lb, ub);
    
    % calculate the posterior values for the pupil fits, given the current
    % measurement and the priors
    pPosteriorTransparent = pPriorSDTransparent.^2.*pInitialFitTransparent./(pPriorSDTransparent.^2+pFitSD.^2) + ...
        pFitSD.^2.*pPriorMeanTransparent./(pPriorSDTransparent.^2+pFitSD.^2);
    
    % re-calculate the fit, fixing the pupil size from the posterior
    lb_pinArea = lb; lb_pinArea(3) = pPosteriorTransparent(3);
    ub_pinArea = ub; ub_pinArea(3) = pPosteriorTransparent(3);
    [pFinalFitTransparent, pFitSD, fitError] = calcPupilLikelihood(Xc,Yc, lb_pinArea, ub_pinArea);
    
    % store results
    pupil.pInitialFitTransparent(i,:) = pInitialFitTransparent';
    pupil.pFinalFitTransparent(i,:) = pFinalFitTransparent';
    pupil.pFitSD(i,:) = pFitSD';
    pupil.fitError(i) = fitError;
    
    % Update the prior, which is the mean of the previous fit
    % values, weighted by a decaying exponential in time
    if i > exponentialTauParam
        range=min([window,i-1]);
        for jj=1:5
            dataVector=squeeze(pupil.pInitialFitTransparent(:,jj))';
            pPriorMeanTransparent(jj) = nansum(exponentialWeights(end-range+1:end).*dataVector(i-range:i-1),2)./nansum(exponentialWeights(end-range+1:end),2);
            pPriorSDTransparent(jj) = nanstd(dataVector(i-range:i-1),exponentialWeights(end-range+1:end));
        end
    end
    
    % plot

    % Plot the pupil boundary data points
    imshow(binP)

    pFitImplicit = ellipse_ex2im(ellipse_transparent2ex(pFinalFitTransparent));
    a = num2str(pFitImplicit(1));
    b = num2str(pFitImplicit(2));
    c = num2str(pFitImplicit(3));
    d = num2str(pFitImplicit(4));
    e = num2str(pFitImplicit(5));
    f = num2str(pFitImplicit(6));
    
    % note that X and Y indices need to be swapped!
    eqt= ['(',a, ')*y^2 + (',b,')*x*y + (',c,')*x^2 + (',d,')*y+ (',e,')*x + (',f,')'];
    
    if isfield(params,'outVideo')
        hold on
        h= ezplot(eqt,[1, 240, 1, 320]);
        set (h, 'Color', 'green')
    end
    
    
    % save frame
    if isfield(params,'outVideo')
        frame   = getframe(ih);
        writeVideo(outObj,frame);
    end
    
end



% save video
if isfield(params,'outVideo')
    close(ih);
    close(outObj);
end
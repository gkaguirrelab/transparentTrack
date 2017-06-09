%% DEV_errorMetric

% looking for a coefficient to assess the accuracy of the tracking.

% CORRELATION strategy: load a *_pupilTrack.mat file. Given the tracking
% params, generate a bw image of the tracked pupil (background white,
% ellipse fitted to the pupil black). Correlate this image with the
% original one. Plot results. Later, do the same with a calibration dot
% video (we expect the correlation to be consistently high in that case).


%% clean

clear all
close all
clc


%% paths
% Set Dropbox directory
%get hostname (for melchior's special dropbox folder settings)
[~,hostname] = system('hostname');
hostname = strtrim(lower(hostname));
if strcmp(hostname,'melchior.uphs.upenn.edu')
    dropboxDir = '/Volumes/Bay_2_data/giulia/Dropbox-Aguirre-Brainard-Lab';
else
    % Get user name
    [~, tmpName] = system('whoami');
    userName = strtrim(tmpName);
    dropboxDir = ['/Users/' userName '/Dropbox-Aguirre-Brainard-Lab'];
end


% for eye tracking
params.projectFolder = 'TOME_data';
params.outputDir = 'TOME_processing';
params.eyeTrackingDir = 'EyeTracking';
params.analysisDir = 'TOME_analysis';

% subject
params.projectSubfolder = 'session2_spatialStimuli';
params.subjectName = 'TOME_3014';
params.sessionDate = '021717';

% run
params.runName = 'tfMRI_RETINO_PA_run01';

ptTrackName = [params.runName '_testTrack.mat'];
videoName = [params.runName '_60hz.avi'];

% calibration
calName = 'GazeCal01_LTcal.mat';


%% load track data
trackData = load (fullfile(dropboxDir,'TOME_processing',params.projectSubfolder,params.subjectName,params.sessionDate,'EyeTracking',ptTrackName));


%% load 60 Hz movie
inVideo = fullfile(dropboxDir,'TOME_processing',params.projectSubfolder,params.subjectName,params.sessionDate,'EyeTracking',videoName);




disp('Loading video file and converting to standard format, may take a couple minutes...');
inObj                   = VideoReader(inVideo);
numFrames               = floor(inObj.Duration*inObj.FrameRate);
grayI                   = zeros([240 320 numFrames],'uint8');

% 
% OVERRIDE NUMFRAMES
numFrames = 3000;

% initialize correlation array
eCorrelation = zeros(1,numFrames);
cCorrelation = zeros(1,numFrames);
pCorrelation = zeros(1,numFrames);

% Convert to gray, resize, crop to livetrack size
for i = 1:numFrames
    thisFrame           = readFrame(inObj);
    tmp                 = rgb2gray(thisFrame);
    tmp2        = imresize(tmp,[486 720]/2); %params.imageSize);
    tmp3 = imcrop(tmp2,[1 1 319 239]);%params.imageCrop);
    grayI(:,:,i) = tmp3;
end
disp('done!');

% ellipse params
ellipseParams = trackData.pupil.ellipseParams;

%% cross correlate every frame with ellipse fit
for ii = 1:numFrames
    
    grayFrame = squeeze(grayI(:,:,ii));
    
    % create a gray image of the ellipse fit
    bwFit = zeros([240 320]);
    
    eParams = ellipseParams(ii);
    if isempty(eParams.status) && ~isempty(eParams.X0_in)
        [Xp, Yp] = calcEllipse(eParams, 360);
        try
            idx = sub2ind([240 320], Xp, Yp);
        catch ME
        end
        if exist ('ME', 'var')
            eCorrelation (ii) = NaN;
            clear ME
            continue
        end
        bwFit(idx) = 1;
        % fill the holes
        bwFit = imfill(bwFit,'holes');
        
        % invert (so that pupil is black)
        bwFit = imcomplement(bwFit);
        
        % convert to gray
        grayFit = uint8(bwFit);
        
        % get correlation value for this frame
        eCorrelation(ii) = (corr2(grayFrame,grayFit));
        
    end
    
end


%% patch correlation
for ii = 1:numFrames
    
    grayFrame = squeeze(grayI(:,:,ii));
    
    % mask the frame with circle params
    % create a mask from circle fitting parameters
    pupilMask = zeros(size(grayFrame));
    pupilMask = insertShape(pupilMask,'FilledCircle',[trackData.pupil.circleX(ii) trackData.pupil.circleY(ii) trackData.pupil.circleRad(ii)*1.5],'Color','white');
    pupilMask = im2bw(pupilMask);
    
    patch = immultiply(grayFrame,pupilMask);
%     patch = imcomplement(patch);
    
    % create a gray image of the ellipse fit
    bwFit = zeros([240 320]);
    
    eParams = ellipseParams(ii);
    if isempty(eParams.status) && ~isempty(eParams.X0_in)
        [Xp, Yp] = calcEllipse(eParams, 360);
        try
            idx = sub2ind([240 320], Xp, Yp);
        catch ME
        end
        if exist ('ME', 'var')
            pCorrelation (ii) = NaN;
            clear ME
            continue
        end
        bwFit(idx) = 1;
        % fill the holes
        bwFit = imfill(bwFit,'holes');
        
        % invert (so that pupil is black)
        bwFit = imcomplement(bwFit);
        
        % mask
        bwFit = immultiply(bwFit,pupilMask);
        
        % convert to gray
        grayFit = uint8(bwFit);
        
        % get correlation value for this frame
        pCorrelation(ii) = (corr2(patch,grayFit));
        
    end
    
end



%% cross correlate with circle fit

for ii = 1:numFrames
    
    grayFrame = squeeze(grayI(:,:,ii));
    
    % create a gray image of the circle fit
    bwFit = zeros([240 320]);
    bwFit = insertShape(bwFit,'FilledCircle',[trackData.pupil.circleX(ii) trackData.pupil.circleY(ii) trackData.pupil.circleRad(ii)*1.5],'Color','white');
    bwFit = im2bw(bwFit);

    % get the fit complement (so that pupil is black
    bwFit = imcomplement(bwFit);
 
    % get correlation value for this frame
    cCorrelation(ii) = (corr2(grayFrame,bwFit));

    
end


%% plot correlation results
fig1 = fullFigure;
plot (eCorrelation)
hold on
plot (cCorrelation)
hold on
plot (trackData.pupil.size ./100)
hold off
xlabel('Frames')
legend ('Ellipse correlation' , 'Circle Fullframe correlation', 'pupil size in px /100 (as tracked by ellipse)')
title([params.subjectName ' ' params.runName ' - FULLFRAME' ] , 'Interpreter', 'none')


% %% patch correlation
% fig2 =  fullFigure;
% plot (pCorrelation)
% hold on
% 
% %% verify that those are blinks
% nanFrames =  find(isnan(eCorrelation));
% 
% for jj = 1: length (nanFrames)
%      grayFrame = squeeze(grayI(:,:,nanFrames(jj)));
%      imshow (grayFrame)
%      pause
% end
% 
% 
% %% verify that those are blinks
% zeroFrames =  find((eCorrelation ==0));
% 
% for jj = 1: length (zeroFrames)
%      grayFrame = squeeze(grayI(:,:,zeroFrames(jj)));
%      imshow (grayFrame)
%      pause   
% end


% Noise/shift characterization using scale calibration videos

% we noticed that the tracked data collected with the liveTrack and our
% customTrack are shifted, even when the tracking is good. This script is
% meant to characterize said shift and give an estimate of the SNR of the
% tracked signals.


%% set paths


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

%% find raw calibration files in data folder
dataDir = fullfile(dropboxDir,'TOME_data',params.projectSubfolder,params.subjectName,params.sessionDate,'EyeTracking');

scaleCalFile = dir(fullfile(dataDir, '*ScaleCal.mat'));

%% find deinterlaced video in processing folder
processingDir = fullfile(dropboxDir,'TOME_processing',params.projectSubfolder,params.subjectName,params.sessionDate,'EyeTracking');

scaleVideos = dir(fullfile(processingDir, '*ScaleCal*60hz.avi'));
%% track deinterlaced video using custom track
for ii = 1 : length(scaleVideos)
    params.acqRate = 60;
    params.inVideo = fullfile(scaleVideos(ii).folder,scaleVideos(ii).name);
    params.outVideo = fullfile(processingDir,[scaleVideos(ii).name '_ellipseTrack.avi']);
    params.outMat = fullfile(processingDir, [scaleVideos(ii).name '_ellipseTrack.mat']);
    params.pupilOnly = 1;
    params.pupilRange   = [30 70];
    params.rangeAdjust  = 0.05;
    params.minMajorAxis = 200; 
    params.maxMajorAxis = 210;
    params.minAspectRatio = 0.8;
     params.rotation = 180; 
     params.rotationSpan = 0.1;
     params.smoothStddev = 1;
     params.randomize = 2;
     params.pupilFit = 'ellipse2';
     params.uniformWeights = false; 
    trackPupil(params);
    
end

%% load tracking data

lt = load (fullfile(dataDir,scaleCalFile.name));
for ii = 1 : length(scaleVideos)
    pt(ii)=load (fullfile(processingDir, [scaleVideos(ii).name '_ellipseTrack.mat']));
end

%% extract liveTrack data

lt5mm= lt.ScaleCal.ReportRaw{3};
lt6mm= lt.ScaleCal.ReportRaw{2};
lt7mm= lt.ScaleCal.ReportRaw{1};

% use all report samples
ct = 0;
for i = 1:length(lt7mm)
    % First field
    ct = ct + 1;
    lt5.pupil.X(ct) = lt5mm(i).PupilCameraX_Ch01;
    lt5.pupil.Y(ct) = lt5mm(i).PupilCameraY_Ch01;
    lt6.pupil.X(ct) = lt6mm(i).PupilCameraX_Ch01;
    lt6.pupil.Y(ct) = lt6mm(i).PupilCameraY_Ch01;
    lt7.pupil.X(ct) = lt7mm(i).PupilCameraX_Ch01;
    lt7.pupil.Y(ct) = lt7mm(i).PupilCameraY_Ch01;
    
    % Second field
    ct = ct + 1;
    lt5.pupil.X(ct) = lt5mm(i).PupilCameraX_Ch02;
    lt5.pupil.Y(ct) = lt5mm(i).PupilCameraY_Ch02;
    lt6.pupil.X(ct) = lt6mm(i).PupilCameraX_Ch02;
    lt6.pupil.Y(ct) = lt6mm(i).PupilCameraY_Ch02;
    lt7.pupil.X(ct) = lt7mm(i).PupilCameraX_Ch02;
    lt7.pupil.Y(ct) = lt7mm(i).PupilCameraY_Ch02;
    
end

%% extract pupil data
pt5.pupil.X = pt(1).pupil.X;
pt5.pupil.Y = pt(1).pupil.Y;
pt6.pupil.X = pt(2).pupil.X;
pt6.pupil.Y = pt(2).pupil.Y;
pt7.pupil.X = pt(3).pupil.X;
pt7.pupil.Y = pt(3).pupil.Y;


%% get median values
% for liveTrack
lt5.median.X = nanmedian(lt5.pupil.X(150:end));
lt5.median.Y = nanmedian(lt5.pupil.Y(150:end));
lt6.median.X = nanmedian(lt6.pupil.X(150:end));
lt6.median.Y = nanmedian(lt6.pupil.Y(150:end));
lt7.median.X = nanmedian(lt7.pupil.X(150:end));
lt7.median.Y = nanmedian(lt7.pupil.Y(150:end));

% for pupilTrack
pt5.median.X = nanmedian(pt5.pupil.X(150:end));
pt5.median.Y = nanmedian(pt5.pupil.Y(150:end));
pt6.median.X = nanmedian(pt6.pupil.X(150:end));
pt6.median.Y = nanmedian(pt6.pupil.Y(150:end));
pt7.median.X = nanmedian(pt7.pupil.X(150:end));
pt7.median.Y = nanmedian(pt7.pupil.Y(150:end));

% get shift
shift5.pupil.X = pt5.median.X - lt5.median.X;
shift6.pupil.X = pt6.median.X - lt6.median.X;
shift7.pupil.X = pt7.median.X - lt7.median.X;

shift5.pupil.Y = pt5.median.Y - lt5.median.Y;
shift6.pupil.Y = pt6.median.Y - lt6.median.Y;
shift7.pupil.Y = pt7.median.Y - lt7.median.Y;

%% plot pupil Y and Y 

% plot
fullFigure;
ylims = [150 200];
xlims = [1 250];
subplot(3,2,1)
plot(lt5.pupil.X);
hold on;
plot(pt5.pupil.X)
grid on
%     ylabel('glint X (normalized)')
%     xlabel('Frames')
legend ('liveTrack','pupilTrack')
title(['pupil X 5 mm']);
% ylim(ylims);
xlim(xlims);
subplot(3,2,2)
plot(lt5.pupil.Y);
hold on;
plot(pt5.pupil.Y)
grid on
%     ylabel('glint X (normalized)')
%     xlabel('Frames')
legend ('liveTrack','pupilTrack')
title(['pupil Y 5 mm']);
% ylim(ylims);
xlim(xlims);

subplot(3,2,3)
plot(lt6.pupil.X);
hold on;
plot(pt6.pupil.X)
grid on
%     ylabel('glint X (normalized)')
%     xlabel('Frames')
legend ('liveTrack','pupilTrack')
title(['pupil X 6 mm']);
% ylim(ylims);
xlim(xlims);
subplot(3,2,4)
plot(lt6.pupil.Y);
hold on;
plot(pt6.pupil.Y)
grid on
%     ylabel('glint X (normalized)')
%     xlabel('Frames')
legend ('liveTrack','pupilTrack')
title(['pupil Y 6 mm']);
% ylim(ylims);
xlim(xlims);

subplot(3,2,5)
plot(lt7.pupil.X);
hold on;
plot(pt7.pupil.X)
grid on
%     ylabel('glint X (normalized)')
%     xlabel('Frames')
legend ('liveTrack','pupilTrack')
title(['pupil X 7 mm']);
% ylim(ylims);
xlim(xlims);
subplot(3,2,6)
plot(lt7.pupil.Y);
hold on;
plot(pt7.pupil.Y)
grid on
%     ylabel('glint X (normalized)')
%     xlabel('Frames')
legend ('liveTrack','pupilTrack')
title(['pupil Y 7 mm']);
% ylim(ylims);
xlim(xlims);



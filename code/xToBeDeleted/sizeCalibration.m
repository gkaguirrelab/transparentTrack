function sizeConversionFactor = sizeCalibration(dropboxDir,params)

%% load files
scaleCalFile = dir(fullfile(dropboxDir, params.projectFolder, params.projectSubfolder, ...
    params.subjectName,params.sessionDate,params.eyeTrackingDir,'*ScaleCal*.mat'));

% load scaleScal file
LT = load(fullfile(scaleCalFile.folder,scaleCalFile.name));

% load all scale Cal Video
rawVids = dir(fullfile(dropboxDir, 'TOME_processing', params.projectSubfolder, ...
    params.subjectName,params.sessionDate,params.eyeTrackingDir,'RawScaleCal*_60hz.avi'));

if isempty(rawVids)
    LTvids = dir(fullfile(dropboxDir, params.projectFolder, params.projectSubfolder, ...
        params.subjectName,params.sessionDate,params.eyeTrackingDir,'*ScaleCal*.avi'));
end


%% track and get conversion factor
if ~isempty(rawVids)
    for rr = 1: length(rawVids)
        fprintf ('\nProcessing calibration %d of %d\n',rr,length(rawVids))
        %get the run name
        params.calName = rawVids(rr).name(1:end-9); %runs
        outDir = fullfile(dropboxDir,'TOME_processing',params.projectSubfolder,params.subjectName,params.sessionDate,'EyeTracking');
        params.acqRate = 60;
        params.pupilFit = 'fixedPupilCut';
        params.ellipseThresh   = [0.94 0.9];
        params.circleThresh = [0.03 0.999];
        params.gammaCorrection = 1;
        params.inVideo = fullfile(outDir,[params.calName '_60hz.avi']);
        params.outVideo = fullfile(outDir,[params.calName '_calTrack.avi']);
        params.outMat = fullfile(outDir, [params.calName '_calTrack.mat']);
        params.pupilOnly = 1;
        params.cutPupil = 1;
        [dotsPX(rr), ~, ~] = trackPupil(params);
    end
    diameters = fliplr(LT.ScaleCal.pupilDiameterMmGroundTruth);
    for rr = 1: length(rawVids)
        PXperMM(rr) = nanmedian(dotsPX(rr).size) / diameters(rr);
    end
    sizeConversionFactor = nanmedian(PXperMM);
else
    % if no raw videos were acquired (early scans), track the livetrack videos
    for rr = 1: length(LTvids)
        fprintf ('\nProcessing calibration %d of %d\n',rr,length(LTvids))
        %get the run name
        params.calName = LTvids(rr).name(1:end-4); %runs
        outDir = fullfile(dropboxDir,'TOME_processing',params.projectSubfolder,params.subjectName,params.sessionDate,'EyeTracking');
        params.acqRate = 10;
        params.pupilFit = 'ellipse';
        params.ellipseThresh   = [0.97 0.9];
        params.circleThresh = [0.05 0.999];
        params.inVideo = fullfile(LTvids(rr).folder, LTvids(rr).name);
        params.outVideo = fullfile(outDir,[params.calName '_calTrack.avi']);
        params.outMat = fullfile(outDir, [params.calName '_calTrack.mat']);
        params.pupilOnly = 1;
        params.cutPupil = 1;
        params.keepOriginalSize = 1;
        [dotsPX(rr), ~, ~] = trackPupil(params);
    end
    diameters = fliplr(LT.ScaleCal.pupilDiameterMmGroundTruth);
    for rr = 1: length(LTvids)
        PXperMM(rr) = nanmedian(dotsPX(rr).size) / diameters(rr);
    end
   sizeConversionFactor = nanmedian(PXperMM);
end % track




%% save out size conversion factor as a mat file
save (fullfile(outDir, 'sizeConversionFactor.mat'), 'sizeConversionFactor')
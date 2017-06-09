function [runParams] = getRunParams (params, dropboxDir)
% header

%% look for scale cal in each date (there is only one Scale calibration file per session)
ScaleCal = dir(fullfile(dropboxDir, params.projectFolder, params.projectSubfolder, ...
    params.subjectName,params.sessionDate,params.eyeTrackingDir,'*ScaleCal*.mat'));

%% look for Gaze Calibration files
GazeCals = dir(fullfile(dropboxDir, params.projectFolder, params.projectSubfolder, ...
    params.subjectName,params.sessionDate,params.eyeTrackingDir,'*LTcal*.mat'));
GazeData = dir(fullfile(dropboxDir, params.projectFolder, params.projectSubfolder, ...
    params.subjectName,params.sessionDate,params.eyeTrackingDir,'*LTdat*.mat'));


%% Create the runParams Struct
runParams.outputDir = params.outputDir;
runParams.projectFolder = params.projectFolder;
runParams.projectSubfolder = params.projectSubfolder;
runParams.subjectName = params.subjectName;
runParams.sessionDate =params.sessionDate;
runParams.eyeTrackingDir = params.eyeTrackingDir;
runParams.runName = params.runName;
runParams.scaleCalName = ScaleCal.name;

%% Assign the appropriate gaze cal file (if gaze files exist)
if ~isempty(GazeCals)
    if length(GazeCals) == 1
        % check if the calibration is good
        LTdatFile = fullfile(dropboxDir, params.projectFolder, params.projectSubfolder, ...
            params.subjectName,params.sessionDate,params.eyeTrackingDir,GazeData.name);
        LTcalFile = fullfile(dropboxDir, params.projectFolder, params.projectSubfolder, ...
            params.subjectName,params.sessionDate,params.eyeTrackingDir,GazeCals.name);
        isGood = checkGazeCal(LTdatFile,LTcalFile);
        % assign gazeCalFile
        if isGood
            runParams.gazeCalName = GazeCals.name;
        else
            warning('No good Gaze Calibration found. Data will not be calibrated for the Gaze')
        end
    else
        % sort Gaze Calibration and Data files by timestamp
        [~,idx] = sort([GazeCals.datenum]);
        [~,idx2] = sort([GazeData.datenum]);
        % check if the calibration are good
        for ii = 1: length(idx)
            LTdatFile = fullfile(dropboxDir, params.projectFolder, params.projectSubfolder, ...
                params.subjectName,params.sessionDate,params.eyeTrackingDir,GazeData(idx2(ii)).name);
            LTcalFile = fullfile(dropboxDir, params.projectFolder, params.projectSubfolder, ...
                params.subjectName,params.sessionDate,params.eyeTrackingDir,GazeCals(idx(ii)).name);
            isGood(idx(ii)) = checkGazeCal(LTdatFile,LTcalFile);
        end
        % get the timestamp for the current run from the corresponding
        % report
        reportFile = dir(fullfile(dropboxDir, params.projectFolder, params.projectSubfolder, ...
            params.subjectName,params.sessionDate,params.eyeTrackingDir,[params.runName '_report.mat']));
        reportTime = reportFile.datenum;
        % take the most recent calibration file acquired before the current run
        for ii = 1: length(idx)
            if GazeCals(idx(ii)).datenum < reportTime && isGood(idx(ii))
                runParams.gazeCalName = GazeCals(idx(ii)).name;
            else
                continue
            end
        end
        % if a preceeding good one is not found, just take the first good
        % one in the session.
        if ~isfield(runParams, 'gazeCalName')
            for ii = 1: length(idx)
                if isGood(idx(ii))
                    runParams.gazeCalName = GazeCals(idx(ii)).name;
                else
                    continue
                end
            end
        end
        % if still no good calibration is found, warn the user
        if ~isfield(runParams, 'gazeCalName')
             warning('No good Gaze Calibration found. Data will not be calibrated for the Gaze')
        end
    end
end


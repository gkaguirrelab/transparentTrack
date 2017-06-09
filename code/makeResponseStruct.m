function [response] = makeResponseStruct(params, dropboxDir)
% this function calibrates and assembles the pupil respose Struct in the standard format.
% 
% Response fields:
%         response.pupilSize 
%         response.gazeEcc 
%         response.gazePol 
%         response.gazeX 
%         response.gazeY  
%         response.timebase 
%         response.metadata

%% get run params
[runParams] = getRunParams (params,dropboxDir);

%% calibrate and asseble response struct
switch params.trackType
    case 'Hybrid'
        
        % calibrate
        calParams.trackType = params.trackType;
        calParams.eyeTrackFile =  fullfile(dropboxDir,runParams.outputDir,runParams.projectSubfolder,...
            runParams.subjectName,runParams.sessionDate,runParams.eyeTrackingDir,...
            [runParams.runName '_pupilTrack.mat']);
        calParams.scaleCalFile = fullfile(dropboxDir,runParams.projectFolder,runParams.projectSubfolder,...
            runParams.subjectName,runParams.sessionDate,runParams.eyeTrackingDir,runParams.scaleCalName);
        % for early session 1 that do not have a Gaze Cal.
        if isfield(runParams,'gazeCalName')
            calParams.gazeCalFile = fullfile(dropboxDir,runParams.projectFolder,runParams.projectSubfolder,...
                runParams.subjectName,runParams.sessionDate,runParams.eyeTrackingDir,runParams.gazeCalName);
        else
            warning('No gaze calibration file found for this session. Will use the first calibration file available from the subject''s session 2.')
            % look for Gaze Calibration files in session 2
            GazeCals = dir(fullfile(dropboxDir, params.projectFolder, params.projectSubfolderTwo, ...
                params.subjectName,params.sessionTwoDate,params.eyeTrackingDir,'*LTcal*.mat'));
            GazeData = dir(fullfile(dropboxDir, params.projectFolder, params.projectSubfolderTwo, ...
                params.subjectName,params.sessionTwoDate,params.eyeTrackingDir,'*LTdat*.mat'));
            % sort Gaze Calibration files by timestamp
            [~,idx] = sort([GazeCals.datenum]);
            [~,idx2] = sort([GazeData.datenum]);
            % check if the calibration are good
            for ii = 1: length(idx)
                LTdatFile = fullfile(dropboxDir, params.projectFolder, params.projectSubfolderTwo, ...
                    params.subjectName,params.sessionTwoDate,params.eyeTrackingDir,GazeData(idx2(ii)).name);
                LTcalFile = fullfile(dropboxDir, params.projectFolder, params.projectSubfolderTwo, ...
                    params.subjectName,params.sessionTwoDate,params.eyeTrackingDir,GazeCals(idx(ii)).name);
                isGood(idx(ii)) = checkGazeCal(LTdatFile,LTcalFile);
            end
            % take the first good GazeCal file
            for ii = 1: length(idx)
                if isGood(idx(ii))
                    calParams.gazeCalFile = fullfile(dropboxDir,runParams.projectFolder,params.projectSubfolderTwo,...
                        runParams.subjectName,params.sessionTwoDate,runParams.eyeTrackingDir,GazeCals(idx(ii)).name);
                    runParams.gazeCalName = GazeCals(idx(ii)).name;
                else
                    continue
                end
            end
        end
        
        % if exists a gaze calibration file by the same name in TOME_processing, load that
        % instead
        calName = runParams.gazeCalName(1:end-10);
        betterCal = dir(fullfile(dropboxDir, 'TOME_processing', params.projectSubfolder, ...
            params.subjectName,params.sessionDate,params.eyeTrackingDir, [calName '_calParams.mat']));
        if ~isempty (betterCal) % && strcmp(params.projectSubfolder, 'session2_spatialStimuli')
            calParams.gazeCalFile = (fullfile (betterCal.folder, betterCal.name));
            runParams.gazeCalName = betterCal.name;
        end
        % if exists a size conversion factor file in TOME_processing, load
        % that instead
        sizeCal = dir(fullfile(dropboxDir, 'TOME_processing', params.projectSubfolder, ...
            params.subjectName,params.sessionDate,params.eyeTrackingDir, 'sizeConversionFactor.mat'));
        if ~isempty (sizeCal)
            calParams.sizeConversionFactor = load(fullfile(sizeCal.folder, sizeCal.name));
        end
            
        
        
        [pupilSize,gaze,pupilError,pupilCut] = calcPupilGaze(calParams);
        
        % load timeBase file
        timeBaseFile = fullfile(dropboxDir,'TOME_processing',runParams.projectSubfolder,...
            runParams.subjectName,runParams.sessionDate,runParams.eyeTrackingDir,...
            [runParams.runName '_timeBase.mat']);
        if exist (timeBaseFile,'file')
            load (timeBaseFile);
        else
            warning ('Run skipped because no timebase file was found');
            response = '';
            return
        end
        
        % assemble response struct values
        response.pupilError = pupilError';
        response.pupilCut = pupilCut';
        response.pupilSize = pupilSize';
        response.gazeEcc = gaze.ecc;
        response.gazePol = gaze.pol;
        response.gazeX = gaze.X';
        response.gazeY = gaze.Y';
        response.timebase = timeBase.pt;
        
        % check that the timebase and the response are the same length
        if length(response.pupilSize)~=length(response.timebase)
            warning ('Timebase and response values are not of the same length')
        end
        
        % metadata
        response.metaData = runParams;
        response.metaData.eyeTrackFile = fullfile(runParams.outputDir,runParams.projectSubfolder,...
            runParams.subjectName,runParams.sessionDate,runParams.eyeTrackingDir,...
            [runParams.runName '_pupilTrack.mat']);
        response.metaData.trackingParams = load(fullfile(dropboxDir,'TOME_processing',runParams.projectSubfolder,...
            runParams.subjectName,runParams.sessionDate,runParams.eyeTrackingDir,...
            [runParams.runName '_trackingParams.mat']));
        response.metaData.scaleCalFile = fullfile(runParams.projectFolder,runParams.projectSubfolder,...
            runParams.subjectName,runParams.sessionDate,runParams.eyeTrackingDir,runParams.scaleCalName);
        response.metaData.gazeCalFile = fullfile(runParams.projectFolder,runParams.projectSubfolder,...
            runParams.subjectName,runParams.sessionDate,runParams.eyeTrackingDir,runParams.gazeCalName);
        response.metaData.trackType = calParams.trackType;
        
        % git info
        % LiveTrack toolbox
        fCheck = which('GetGitInfo');
        if ~isempty(fCheck)
            thePath = fileparts(mfilename('fullpath'));
            gitInfo = GetGITInfo(thePath);
        else
            gitInfo = 'function ''GetGITInfo'' not found';
        end
        response.metaData.gitInfo = gitInfo;
        
end

function saveGazePostData( dataRootDir, dataSaveDir, varargin )
% Saves plots of the quality of ellipse fitting to the pupil perimeter
%
% Syntax:
%  saveGazePostData( dataRootDir, dataSaveDir )
%
% Description:
%
%
% Inputs:
%   dataRootDir           - Full path to a directory that contains pupil
%                           data files, and/or subdirectories containing
%                           these files.
%   dataSaveDir           - Where to save the gaze data.
%
% Optional key/value pairs:
%
% Outputs:
%   none
%
% Examples:
%{
    dropboxBaseDir=fullfile(getpref('eyeTrackTOMEAnalysis','dropboxBaseDir'));
    dataRootDir=fullfile(dropboxBaseDir,'TOME_processing','session2_spatialStimuli');
    dataSaveDir=fullfile(dataRootDir,'pupilDataQAPlots_eyePose_MOVIE_July2020');
    saveGazePostData( dataRootDir, dataSaveDir,'acquisitionStem','tfMRI_MOVIE')
%}
%{
    dropboxBaseDir=fullfile(getpref('eyeTrackTOMEAnalysis','dropboxBaseDir'));
    dataRootDir=fullfile(dropboxBaseDir,'TOME_processing','session2_spatialStimuli');
    dataSaveDir=fullfile(dataRootDir,'pupilDataQAPlots_eyePose_RETINO_July2020');
    saveGazePostData( dataRootDir, dataSaveDir,'acquisitionStem','tfMRI_RETINO')
%}

%% input parser
p = inputParser;

% Required
p.addRequired('dataRootDir',@ischar);
p.addOptional('dataSaveDir',[],@(x)(isempty(x) || ischar(x)));

% Optional
p.addParameter('rmseThreshold',3,@isscalar);
p.addParameter('acquisitionStem','rfMRI_REST',@ischar);


% parse
p.parse(dataRootDir,dataSaveDir,varargin{:})


% Data window
deltaT = 1000/60;
dataWindowDurMsec = 336 * 1000 - deltaT;

% Obtain the paths to all of the pupil data files within the specified
% directory, including within sub-drectories.
fileListStruct=dir(fullfile(dataRootDir,'**',[p.Results.acquisitionStem '*_pupil.mat']));

% If we found at least one pupil data file, then proceed.
if ~isempty(fileListStruct)
    
    % Get a list of the directories that contain a pupil data file. We will
    % treat each of these as a session, and the pupil data files within as
    % different acquisitions from that session.
    fileListCell=struct2cell(fileListStruct);
    uniqueDirNames=unique(fileListCell(2,:));
    
    % Create a list of name tags that are the unique character components
    % of the uniqueDirNames. We trim off all characters that are shared
    % across the full path of the dir names.
    if length(uniqueDirNames) > 1
        nameTags = uniqueDirNames;
        stillTrimming = true;
        while stillTrimming
            thisChar=cell2mat(cellfun(@(x) double(x(1)),nameTags,'UniformOutput',false));
            if range(thisChar)==0
                nameTags=cellfun(@(x) x(2:end),nameTags,'UniformOutput',false);
            else
                stillTrimming=false;
            end
        end
        stillTrimming = true;
        while stillTrimming
            thisChar=cell2mat(cellfun(@(x) double(x(end)),nameTags,'UniformOutput',false));
            if range(thisChar)==0
                nameTags=cellfun(@(x) x(1:end-1),nameTags,'UniformOutput',false);
            else
                stillTrimming=false;
            end
        end
    end
    
    % Replace file system delimeters with valid characters. That way, we
    % can use these name tags as filenames for the resulting plots.
    nameTags=cellfun(@(x) strrep(x,filesep,'_'),nameTags,'UniformOutput',false);
    
    gazeData = {};
    
    % Loop through the set of sessions
    for ii = 1:length(uniqueDirNames)
        
        % Get the list of acquisition file names
        acqList = find(strcmp(fileListCell(2,:),uniqueDirNames{ii}));
        
        % Loop through the acquisitions
        for jj = 1: length(acqList)
            
            % Obtain the path to this pupil data file
            pupilFullFileName = fullfile(fileListCell{2,acqList(jj)},fileListCell{1,acqList(jj)});
            load(pupilFullFileName,'pupilData');
            
            % Grab just the filename for this pupil data, omitting the
            % path. We will use this to label the plot
            [pupilFilePath, pupilFileName] = fileparts(pupilFullFileName);
            fileNameStem = strsplit(pupilFileName,'_pupil');
            fileNameStem = fileNameStem{1};
            
            % Load the glint
            glintbaseFileName = fullfile(pupilFilePath,[fileNameStem,'_glint.mat']);
            load(glintbaseFileName,'glintData');
            
            % Load the sceneGeometry
            sceneGeomFileName = fullfile(pupilFilePath,[fileNameStem,'_sceneGeometry.mat']);
            if ~isfile(sceneGeomFileName)
                continue
            end
            load(sceneGeomFileName,'sceneGeometry');
            
            % If the sceneGeometry lacks a pose registration, skip forward
            if ~isfield(sceneGeometry.screenPosition.poseRegParams,'R')
                sceneGeomFileName
                continue
            end
            
            % Get and store any spectacle magnification for this scene
            spectacleMag = 1;
            if isfield(sceneGeometry.refraction.retinaToCamera.magnification,'spectacle')
                spectacleMag = sceneGeometry.refraction.retinaToCamera.magnification.spectacle;
            end            
            
            % Find the frames that lack a glint
            noGlint = true(size(pupilData.initial.ellipses.RMSE));
            stopFrame = min([length(glintData.X),length(pupilData.initial.ellipses.RMSE)]);
            noGlint(1:stopFrame) = isnan(glintData.X(1:stopFrame));
            
            % Load the timebase file
            timebaseFileName = fullfile(pupilFilePath,[fileNameStem,'_timebase.mat']);
            load(timebaseFileName,'timebase');
            
            % Check that there is a radiusSmoothed field; otherwise
            % continue
            if ~isfield(pupilData,'radiusSmoothed')
                continue
            end
            
            % Obtain the vector of good and bad time points for the
            % radius smoothed
            highRMSE = pupilData.radiusSmoothed.ellipses.RMSE > p.Results.rmseThreshold;
            nanFit = isnan(pupilData.radiusSmoothed.ellipses.RMSE);
            fitAtBound = pupilData.radiusSmoothed.eyePoses.fitAtBound;
            goodRadiusSmoothed = logical(~highRMSE(1:length(fitAtBound)) .* ~fitAtBound(1:length(fitAtBound)) .* ~noGlint(1:length(fitAtBound)) .* ~nanFit(1:length(fitAtBound)));
            
            % If there are fewer than 66% good points, skip this
            % acquisition
            if sum(goodRadiusSmoothed)<(length(goodRadiusSmoothed)/1.5)
                continue
            end
            
            % Convert the eyePose to gaze position
            f = sceneGeometry.screenPosition.poseRegParams.R * [pupilData.radiusSmoothed.eyePoses.values(:,1), pupilData.radiusSmoothed.eyePoses.values(:,2)]' + sceneGeometry.screenPosition.poseRegParams.t;
            
            % Account for the effect of spectacle magnification
            f = f./spectacleMag;
            
            % Add the pupil
            f(3,:) = pupilData.radiusSmoothed.eyePoses.values(:,4);
            
            % Nan out the bad points
            f(:,~goodRadiusSmoothed)=nan;
            
            % Retain just the portion that is in the time frame, and
            % resample
            [~, startFrame] = min(abs(timebase.values));
            [~, endFrame] = min(abs(timebase.values-dataWindowDurMsec));
            
            if size(f,2)<endFrame
                f(:,end:endFrame)=nan;
            end
            
            for vv = 1:3
               vq(vv,:) = interp1(timebase.values(startFrame:endFrame),f(vv,startFrame:endFrame),0:deltaT:dataWindowDurMsec);
            end
            
            rmse = interp1(timebase.values(startFrame:endFrame),pupilData.radiusSmoothed.ellipses.RMSE(startFrame:endFrame),0:deltaT:dataWindowDurMsec);
            
            if ~isfield(gazeData,fileNameStem)
                gazeData.timebase = 0:deltaT:dataWindowDurMsec;
                gazeData.(fileNameStem).vq(1,:,:) = vq;
                gazeData.(fileNameStem).RMSE(1,:) = rmse;
                gazeData.(fileNameStem).spectacleMag(1) = spectacleMag;
                gazeData.(fileNameStem).nameTags{1} = nameTags{ii};
                gazeData.(fileNameStem).filePathStem{1} = fullfile(pupilFilePath,fileNameStem);
            else
                gazeData.(fileNameStem).vq(end+1,:,:) = vq;
                gazeData.(fileNameStem).RMSE(end+1,:) = rmse;
                gazeData.(fileNameStem).spectacleMag(end+1) = spectacleMag;
                gazeData.(fileNameStem).nameTags{end+1} = nameTags{ii};
                gazeData.(fileNameStem).filePathStem{end+1} = fullfile(pupilFilePath,fileNameStem);
            end
            
        end % loop over acquisitions       
                
    end % loop over sessions
    
    dataFileName = fullfile(dataSaveDir,'gazeData.mat');
    save(dataFileName,'gazeData');
    

end % we have at least one session

end % saveGazePostData


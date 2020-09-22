function plotRelCameraPosition( dataRootDir, plotSaveDir, varargin )
% Saves plots of the quality of ellipse fitting to the pupil perimeter
%
% Syntax:
%  plotRelCameraPosition( dataRootDir, plotSaveDir )
%
% Description:
%   Given a directory, the routine will perform a recursive search through
%   sub-directories to find all files with the suffix "_pupil.mat".
%   Multiple pupil data files within a sub-directory will be treated as
%   different acquisitions from a single session. A plot is generated for
%   each session that presents measures of camera translation for each
%   acquisition, for each specified stage. Each session results in a
%   separate plot, all of which are saved in the directory specified by
%   "plotSaveDir". If the plotSaveDir variable is empty, the plots will be
%   displayed and not saved.
%
%   Interpretation of the plot elements:
%     - Time is given in minutes relative to the start of the fMRI scan
%     - Gray, red, and blue are translation in horizontal, vertical, and
%     depth.
%
% Inputs:
%   dataRootDir           - Full path to a directory that contains pupil
%                           data files, and/or subdirectories containing
%                           these files.
%   plotSaveDir           - Optional. Full path the a directory within
%                           which the QA plots should be saved. If not
%                           specified, or passed as empty, the plots will
%                           be displayed but not saved.
%
% Optional key/value pairs:
%
% Outputs:
%   none
%
% Examples:
%{
    dropboxBaseDir=fullfile(getpref('eyeTrackTOMEAnalysis','dropboxBaseDir'));
    dataRootDir=fullfile(dropboxBaseDir,'TOME_processing','session1_restAndStructure');
    dataSaveDir=fullfile(dataRootDir,'pupilDataQAPlots_cameraTrans_REST_July2020');
    plotRelCameraPosition( dataRootDir,dataSaveDir,'acquisitionStem','rfMRI_REST','nColumns',4,'selectSubjects',{'3001'});
%}
%{
    dropboxBaseDir=fullfile(getpref('eyeTrackTOMEAnalysis','dropboxBaseDir'));
    dataRootDir=fullfile(dropboxBaseDir,'TOME_processing','session2_spatialStimuli');
    dataSaveDir=fullfile(dataRootDir,'pupilDataQAPlots_cameraTrans_FLASH_July2020');
    plotRelCameraPosition( dataRootDir, dataSaveDir,'acquisitionStem','tfMRI_FLASH','nColumns',2)
%}
%{
    dropboxBaseDir=fullfile(getpref('eyeTrackTOMEAnalysis','dropboxBaseDir'));
    dataRootDir=fullfile(dropboxBaseDir,'TOME_processing','session2_spatialStimuli');
    dataSaveDir=fullfile(dataRootDir,'pupilDataQAPlots_cameraTrans_RETINO_July2020');
    plotRelCameraPosition( dataRootDir, dataSaveDir,'acquisitionStem','tfMRI_RETINO','nColumns',4)
%}
%{
    dropboxBaseDir=fullfile(getpref('eyeTrackTOMEAnalysis','dropboxBaseDir'));
    dataRootDir=fullfile(dropboxBaseDir,'TOME_processing','session2_spatialStimuli');
    dataSaveDir=fullfile(dataRootDir,'pupilDataQAPlots_cameraTrans_MOVIE_July2020');
    plotRelCameraPosition( dataRootDir, dataSaveDir,'acquisitionStem','tfMRI_MOVIE','nColumns',4,'selectSubjects',{'3022'})
%}
%{
    % Special for the subjects without gazeCalibration data
    dropboxBaseDir=fullfile(getpref('eyeTrackTOMEAnalysis','dropboxBaseDir'));
    dataRootDir=fullfile(dropboxBaseDir,'TOME_processing','session2_spatialStimuli');
    dataSaveDir=fullfile(dataRootDir,'pupilDataQAPlots_cameraTrans_RETINO_July2020');
    plotRelCameraPosition( dataRootDir, dataSaveDir,'acquisitionStem','tfMRI_RETINO','nColumns',4,'selectSubjects',{'3001','3002','3003','3005'})
%}
%{
    % Special for the subjects without gazeCalibration data
    dropboxBaseDir=fullfile(getpref('eyeTrackTOMEAnalysis','dropboxBaseDir'));
    dataRootDir=fullfile(dropboxBaseDir,'TOME_processing','session2_spatialStimuli');
    dataSaveDir=fullfile(dataRootDir,'pupilDataQAPlots_cameraTrans_MOVIE_July2020');
    plotRelCameraPosition( dataRootDir, dataSaveDir,'acquisitionStem','tfMRI_MOVIE','nColumns',4,'stagesToPlot',{'initial','radiusSmoothed'},'selectSubjects',{'3001','3002','3003','3005'})
%}

%% input parser
p = inputParser;

% Required
p.addRequired('dataRootDir',@ischar);
p.addOptional('plotSaveDir',[],@(x)(isempty(x) || ischar(x)));

% Optional
p.addParameter('rmseThreshold',1.5,@isscalar);
p.addParameter('stagesToPlot',{'initial','estimateSceneParams','radiusSmoothed'},@iscell);
p.addParameter('transDirColors',{[0.25 0.25 0.25],[1 0.25 0.25],[0.25 0.25 1]},@iscell);
p.addParameter('yRangeIncrement',[5 5 0.25],@isnumeric);
p.addParameter('xLim',[-0.5 5.6],@isnumeric);
p.addParameter('yAxisLabels',{'initial trans [mm]','scene estimation [mm]','final trans [mm]'},@iscell);
p.addParameter('nColumns',4,@isscalar);
p.addParameter('acquisitionStem','rfMRI_REST',@ischar);
p.addParameter('selectSubjects',{},@iscell);


% parse
p.parse(dataRootDir,plotSaveDir,varargin{:})

% A constant that we will use later
msecToMin = 1/1000/60;

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
    
    % Loop through the set of sessions
    for ii = 1:length(uniqueDirNames)
        
        % If selectSubjects is defined, check if this name is on the list
        if ~isempty(p.Results.selectSubjects)
            if ~contains(uniqueDirNames{ii},p.Results.selectSubjects)
                continue
            end
        end
        
        % Get the list of acquisition file names
        acqList = find(strcmp(fileListCell(2,:),uniqueDirNames{ii}));
        
        % Prepare a figure. If a plotSaveDir was not specified, then
        % display the figures. Otherwise, hide them.
        if isempty(plotSaveDir)
            figHandle=figure('visible','on','Name',nameTags{ii});
        else
            figHandle=figure('visible','on','Name',nameTags{ii});
        end
        
        % Format the figure
        set(gcf,'PaperOrientation','landscape');
        set(figHandle, 'Units','inches')
        height = 6;
        width = 11;
        set(figHandle, 'Position',[25 5 width height],...
            'PaperSize',[width height],...
            'PaperPositionMode','auto',...
            'renderer','painters', ...            '
            'Color','w' ...
            );
        
        
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
            
            % Load the camera position
            relCamPosFullFileName = strrep(fullfile(fileListCell{2,acqList(jj)},fileListCell{1,acqList(jj)}),'_pupil.mat','_relativeCameraPosition.mat');
            load(relCamPosFullFileName,'relativeCameraPosition');
                        
            % Load the glint
            timebaseFileName = fullfile(pupilFilePath,[fileNameStem,'_glint.mat']);
            load(timebaseFileName,'glintData');
            
            % Find the frames that lack a glint
            noGlint = true(size(pupilData.initial.ellipses.RMSE));
            stopFrame = min([length(glintData.X),length(pupilData.initial.ellipses.RMSE)]);
            noGlint(1:stopFrame) = isnan(glintData.X(1:stopFrame));
            
            % Load the timebase file
            timebaseFileName = fullfile(pupilFilePath,[fileNameStem,'_timebase.mat']);
            load(timebaseFileName,'timebase');
            
            % Check that the last desired field is present, otherwise
            % continue
            lastField = p.Results.stagesToPlot{end};
            if ~isfield(relativeCameraPosition,lastField)
                continue
            end

            % Obtain the vector of good and bad time points for last field
            if isfield(pupilData,lastField)
                highRMSE = pupilData.(lastField).ellipses.RMSE > p.Results.rmseThreshold;
                fitAtBound = pupilData.(lastField).eyePoses.fitAtBound;
                goodIdx = logical(~highRMSE .* ~fitAtBound .* ~noGlint(1:length(fitAtBound)));
            end
            
            % Loop over the fields to plot
            for kk=1:length(p.Results.stagesToPlot)
                
                % Define the subplot for this acqusition
                subplot(length(p.Results.stagesToPlot),p.Results.nColumns,(kk-1)*p.Results.nColumns+jj,'align');
                
                % Plot the three translation directions
                for dd = 1:3
                    vec = relativeCameraPosition.(p.Results.stagesToPlot{kk}).values(dd,:);
                    if strcmp(p.Results.stagesToPlot{kk},'radiusSmoothed')
                        vec(~goodIdx)=nan;
                    end
                    plot(timebase.values(1:length(fitAtBound))*msecToMin,vec(1:length(fitAtBound)),'Color',p.Results.transDirColors{dd},'LineWidth',0.25);
                    hold on
                end
                
                % Set the plot limits
                xlim(p.Results.xLim);
                xticks(fix(p.Results.xLim(1)):1:fix(p.Results.xLim(2)))
                ylim([-5 5]);
                
                % Remove the chart junk
                if jj ~= 1
                    set(gca,'YColor','none','TickDir','out')
                else
                    % Add a y-axis label
                    ylabel(p.Results.yAxisLabels{kk});
                end
                set(gca,'TickDir','out')
                if kk == 1
                    if isfield(pupilData,'radiusSmoothed')
                        title({fileNameStem,pupilData.radiusSmoothed.meta.timestamp},'Interpreter', 'none');
                    else
                        if isfield(pupilData,'radiusSmoothed')
                            title({fileNameStem,pupilData.sceneConstrained.meta.timestamp},'Interpreter', 'none');
                        else
                            title(fileNameStem,'Interpreter', 'none');
                        end
                    end
                end
                if kk ~= length(p.Results.stagesToPlot)
                    set(gca,'XColor','none')
                else
                    if jj==1
                        xlabel('time from scan start [mins]');
                    end
                end
                
                box off
                
            end % loop over eyePose params
        end % loop over acquisitions
        
        % Save the plot if a plotSaveDir has been defined
        if ~isempty(plotSaveDir)
            plotFileName =fullfile(plotSaveDir,[nameTags{ii} '_cameraTrans.png']);
            print(gcf,plotFileName,'-dpng','-r600')
            close(figHandle)
            
        end
        
    end % loop over sessions
end % we have at least one session

end % plotPupilDataEyePose

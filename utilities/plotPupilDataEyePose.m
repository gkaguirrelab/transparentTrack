function plotPupilDataEyePose( dataRootDir, plotSaveDir, varargin )
% Saves plots of the quality of ellipse fitting to the pupil perimeter
%
% Syntax:
%  plotPupilDataEyePose( dataRootDir, plotSaveDir )
%
% Description:
%   Given a directory, the routine will perform a recursive search through
%   sub-directories to find all files with the suffix "_pupil.mat".
%   Multiple pupil data files within a sub-directory will be treated as
%   different acquisitions from a single session. A plot is generated for
%   each session that presents measures of the quality of the initial
%   ellipse fit to the pupil perimeter for each acquisition. Each session
%   results in a separate plot, all of which are saved in the directory
%   specified by "plotSaveDir". If the plotSaveDir variable is empty, the
%   plots will be displayed and not saved.
%
%   Interpretation of the plot elements:
%     - Time is given in minutes relative to the start of the fMRI scan
%     - The gray line is the eyePose fit value for the sceneConstrained
%       pupilData
%     - The red plot points are the sceneConstrained eyePose fit values,
%       but only for those points that meet criteria for having a "good"
%       fit.
%     - The black line is the radiusSmoothed eyePose fit values.
%     - For a given frame to be considered to have a "good" fit, the
%       following criteria must be met:
%         - the frame must have a glint
%         - the RMSE of the fit to the pupil perimeter points (in units of
%           pixels) must be below the rmseThreshold value (default = 1.5)
%         - the fit for the frame must not have hit the upper or lower
%           bounds on the eyePose params
%     - At the bottom of the plot, gray dots mark frames with high RMSE,
%       and blue dots mark points where the fit hit the eyePose bounds.
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
%  'histMax'              - The upper bound on the last bin that is used to
%                           construct the histogram of RMSE values.
%  'numPlotRows'          - The number of acquisition rows that are
%                           displayed on a plot. This should be set equal
%                           to the maximum number of acquisitions that are
%                           expected across all sessions.
%  'fitLabel'             - The field of pupilData for which the RMSE
%                           ellipse fit values are to be analyzed. The
%                           default value is for the initial ellipse fit.
%
% Outputs:
%   none
%
% Examples:
%{
    dropboxBaseDir=fullfile(getpref('eyeTrackTOMEAnalysis','dropboxBaseDir'));
    dataRootDir=fullfile(dropboxBaseDir,'TOME_processing','session2_spatialStimuli');
    dataSaveDir=fullfile(dataRootDir,'pupilDataQAPlots_eyePose_FLASH_July2020');
    plotPupilDataEyePose( dataRootDir, dataSaveDir,'acquisitionStem','tfMRI_FLASH','nColumns',2,'selectSubjects',{'TOME_3025'})
%}
%{
    dropboxBaseDir=fullfile(getpref('eyeTrackTOMEAnalysis','dropboxBaseDir'));
    dataRootDir=fullfile(dropboxBaseDir,'TOME_processing','session2_spatialStimuli');
    dataSaveDir=fullfile(dataRootDir,'pupilDataQAPlots_eyePose_MOVIE_July2020');
    plotPupilDataEyePose( dataRootDir, dataSaveDir,'acquisitionStem','tfMRI_MOVIE','nColumns',4,'selectSubjects',{'TOME_3022'})
%}
%{
    dropboxBaseDir=fullfile(getpref('eyeTrackTOMEAnalysis','dropboxBaseDir'));
    dataRootDir=fullfile(dropboxBaseDir,'TOME_processing','session2_spatialStimuli');
    dataSaveDir=fullfile(dataRootDir,'pupilDataQAPlots_eyePose_RETINO_July2020');
    plotPupilDataEyePose( dataRootDir, dataSaveDir,'acquisitionStem','tfMRI_RETINO','nColumns',4,'selectSubjects',{'TOME_3011'})
%}
%{
    dropboxBaseDir=fullfile(getpref('eyeTrackTOMEAnalysis','dropboxBaseDir'));
    dataRootDir=fullfile(dropboxBaseDir,'TOME_processing','session1_restAndStructure');
    dataSaveDir=fullfile(dataRootDir,'pupilDataQAPlots_eyePose_July2020');
    plotPupilDataEyePose( dataRootDir, dataSaveDir,'acquisitionStem','rfMRI_REST','nColumns',4,'rmseThreshold',2.25,'selectSubjects',{'TOME_3003'})
%}


%% input parser
p = inputParser;

% Required
p.addRequired('dataRootDir',@ischar);
p.addOptional('plotSaveDir',[],@(x)(isempty(x) || ischar(x)));

% Optional
p.addParameter('rmseThreshold',1.5,@isscalar);
p.addParameter('eyePoseParamsToPlot',[1 2 4],@isnumeric);
p.addParameter('yRangeIncrement',[5 5 0.25],@isnumeric);
p.addParameter('xLim',[-0.5 5.6],@isnumeric);
p.addParameter('yAxisLabels',{'azimuth [deg]','elevation [deg]','radius [mm]'},@iscell);
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
        
        % First loop through the acquisitions to determine the y-axis range
        % for each eye parameter.
        lb = []; ub = [];
        for mm = 1:length(p.Results.eyePoseParamsToPlot)
            lb_byParam = nan(1,length(p.Results.eyePoseParamsToPlot));
            ub_byParam = nan(1,length(p.Results.eyePoseParamsToPlot));
            for nn = 1:length(acqList)
                % Obtain the path to this pupil data file and load it
                pupilFullFileName = fullfile(fileListCell{2,acqList(nn)},fileListCell{1,acqList(nn)});
                load(pupilFullFileName,'pupilData');
                glintFullFileName = strrep(fullfile(fileListCell{2,acqList(nn)},fileListCell{1,acqList(nn)}),'_pupil.mat','_glint.mat');
                load(glintFullFileName,'glintData');
                
                % Find the frames that lack a glint
                noGlint = true(size(pupilData.initial.ellipses.RMSE));
                stopFrame = min([length(glintData.X),length(pupilData.initial.ellipses.RMSE)]);
                noGlint(1:stopFrame) = isnan(glintData.X(1:stopFrame));
                
                % Bounds are set for the sceneConstrained
                if ~isfield(pupilData,'sceneConstrained')
                    continue
                else
                    highRMSE = pupilData.sceneConstrained.ellipses.RMSE > p.Results.rmseThreshold;
                    fitAtBound = pupilData.sceneConstrained.eyePoses.fitAtBound;
                    goodSceneConstrained = logical(~highRMSE .* ~fitAtBound .* ~noGlint(1:length(fitAtBound)));
                    
                    % If there are no good data points, fill withn ans
                    if sum(goodSceneConstrained)==0
                        lb_byParam(nn) = nan;
                        ub_byParam = nan;
                    else
                        lb_byParam(nn) = floor(min(pupilData.sceneConstrained.eyePoses.values(goodSceneConstrained,p.Results.eyePoseParamsToPlot(mm))) ./ p.Results.yRangeIncrement(mm)).*p.Results.yRangeIncrement(mm);
                        ub_byParam(nn) = ceil(max(pupilData.sceneConstrained.eyePoses.values(goodSceneConstrained,p.Results.eyePoseParamsToPlot(mm))) ./ p.Results.yRangeIncrement(mm)).*p.Results.yRangeIncrement(mm);
                    end
                end
            end
            lb(mm) = nanmin(lb_byParam);
            ub(mm) = nanmax(ub_byParam);
        end
        
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
            timebaseFileName = fullfile(pupilFilePath,[fileNameStem,'_glint.mat']);
            load(timebaseFileName,'glintData');
            
            % Find the frames that lack a glint
            noGlint = true(size(pupilData.initial.ellipses.RMSE));
            stopFrame = min([length(glintData.X),length(pupilData.initial.ellipses.RMSE)]);
            noGlint(1:stopFrame) = isnan(glintData.X(1:stopFrame));
            
            % Load the timebase file
            timebaseFileName = fullfile(pupilFilePath,[fileNameStem,'_timebase.mat']);
            load(timebaseFileName,'timebase');
            
            % Check that there is at least a sceneConstrained field;
            % otherwise continue
            if ~isfield(pupilData,'sceneConstrained')
                continue
            end
            
            % Obtain the vector of good and bad time points for the
            % sceneConstrained
            highRMSE = pupilData.sceneConstrained.ellipses.RMSE > p.Results.rmseThreshold;
            fitAtBound = pupilData.sceneConstrained.eyePoses.fitAtBound;
            goodSceneConstrained = logical(~highRMSE .* ~fitAtBound .* ~noGlint(1:length(fitAtBound)));
            
            % Obtain the vector of good and bad time points for the
            % radius smoothed, if that exists
            if isfield(pupilData,'radiusSmoothed')
                highRMSE = pupilData.radiusSmoothed.ellipses.RMSE > p.Results.rmseThreshold;
                fitAtBound = pupilData.radiusSmoothed.eyePoses.fitAtBound;
                goodRadiusSmoothed = logical(~highRMSE .* ~fitAtBound .* ~noGlint(1:length(fitAtBound)));
            end
            
            % Loop over the 3 eyePose parameters to be plotted
            for kk=1:length(p.Results.eyePoseParamsToPlot)
                
                % Define the subplot for this acqusition
                subplot(length(p.Results.eyePoseParamsToPlot),p.Results.nColumns,(kk-1)*p.Results.nColumns+jj,'align');
                
                % Plot all values from the sceneConstrained time-series as
                % a gray line
                plot(timebase.values(1:length(goodRadiusSmoothed))*msecToMin,pupilData.sceneConstrained.eyePoses.values(:,p.Results.eyePoseParamsToPlot(kk)),'-','Color',[0.95 0.95 0.95],'LineWidth',0.5);
                hold on
                
                % Now just the "good" sceneConstrained values, plotted as
                % transparent red points.
                if sum(goodSceneConstrained)>0
                    hLineRed = plot(timebase.values(goodSceneConstrained)*msecToMin,pupilData.sceneConstrained.eyePoses.values(goodSceneConstrained,p.Results.eyePoseParamsToPlot(kk)),'o','MarkerSize',1);
                    drawnow
                    hMarkerRed = hLineRed.MarkerHandle;
                    hMarkerRed.FaceColorData = uint8(255*[1; 0; 0; 0.025]);
                    hMarkerRed.FaceColorType = 'truecoloralpha';
                    hMarkerRed.EdgeColorData = uint8([0; 0; 0; 0]);
                    
                    % Plot the radiusSmoothed time-series as a thin black line
                    if isfield(pupilData,'radiusSmoothed')
                        vec = pupilData.radiusSmoothed.eyePoses.values(:,p.Results.eyePoseParamsToPlot(kk));
                        vec(~goodRadiusSmoothed)=nan;
                        hLineBlack = plot(timebase.values(1:length(vec))*msecToMin,vec,'-k','LineWidth',0.25);
                        hLineBlack.Color(4) = 0.5;
                    end
                    
                    % Add markers for high RMSE plot points
                    lowY = lb(kk) + (ub(kk)-lb(kk))/20;
                    yPosition = repmat(lowY,size(timebase.values(highRMSE)));
                    hLineGray = plot(timebase.values(highRMSE)/1000/60,yPosition,'o','MarkerSize',0.75);
                    drawnow
                    if ~isempty(hLineGray)
                        hMarkerGray = hLineGray.MarkerHandle;
                        hMarkerGray.FaceColorData = uint8(255*[0; 0; 0; 0.75]);
                        hMarkerGray.FaceColorType = 'truecoloralpha';
                        hMarkerGray.EdgeColorData = uint8([0; 0; 0; 0]);
                    end
                    
                    % Add markers for at bound plot points
                    lowY = lb(kk) + (ub(kk)-lb(kk))/15;
                    yPosition = repmat(lowY,size(timebase.values(fitAtBound)));
                    if isfield(pupilData.sceneConstrained.eyePoses,'fitAtBound')
                        hLineBlue = plot(timebase.values(fitAtBound)/1000/60,yPosition,'o','MarkerSize',0.75);
                        drawnow
                        if ~isempty(hLineBlue)
                            hMarkerBlue = hLineBlue.MarkerHandle;
                            hMarkerBlue.FaceColorData = uint8(255*[0; 0; 0.75; 1]);
                            hMarkerBlue.FaceColorType = 'truecoloralpha';
                            hMarkerBlue.EdgeColorData = uint8([0; 0; 0; 0]);
                        end
                    end
                    
                    % Set the plot limits
                    xlim(p.Results.xLim);
                    xticks(fix(p.Results.xLim(1)):1:fix(p.Results.xLim(2)))
                    ylim([lb(kk) ub(kk)]);
                end
                
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
                        title({fileNameStem,pupilData.sceneConstrained.meta.timestamp},'Interpreter', 'none');
                    end
                end
                if kk ~= length(p.Results.eyePoseParamsToPlot)
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
            plotFileName =fullfile(plotSaveDir,[nameTags{ii} '_eyePose.png']);
            print(gcf,plotFileName,'-dpng','-r600')
            close(figHandle)
            
        end
        
    end % loop over sessions
end % we have at least one session

end % plotPupilDataEyePose

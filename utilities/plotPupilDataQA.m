function plotPupilDataQA( dataRootDir, plotSaveDir, varargin )



%% input parser
p = inputParser;

% Required
p.addRequired('dataRootDir',@ischar);
p.addRequired('plotSaveDir',@ischar);

% Optional
p.addParameter('histMax',2,@isnumeric);
p.addParameter('numPlotRows',16,@isnumeric);
p.addParameter('fitLabel','initial',@ischar);

% parse
p.parse(dataRootDir,plotSaveDir,varargin{:})

histEdges = 0:0.05:p.Results.histMax;

fileListStruct=subdir(fullfile(dataRootDir,'*_pupil.mat'));

if ~isempty(fileListStruct)
    fileListCell=struct2cell(fileListStruct);
    uniqueDirNames=unique(fileListCell(2,:));
    % Create a list of name tags that are the unique character components
    % of the uniqueDirNames
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
    % Replace file system delimeters with valid characters
    nameTags=cellfun(@(x) strrep(x,filesep,'_'),nameTags,'UniformOutput',false);
    
    for ii = 1:length(uniqueDirNames)
        acqList = find(strcmp(fileListCell(2,:),uniqueDirNames{ii}));
        figHandle=figure('visible','off');
        plotPos = [];
        set(gcf,'PaperOrientation','landscape');
        set(figHandle, 'Units','inches')
        height = 6;
        width = 11;
        for jj = 1: length(acqList)
            pupilFullFileName = fullfile(dataRootDir,fileListCell{1,acqList(jj)});
            dataLoad=load(pupilFullFileName);
            pupilData=dataLoad.pupilData;
            clear dataLoad
            if isfield(pupilData,p.Results.fitLabel)
                sp=subplot(p.Results.numPlotRows+1,1,jj,'align');
                errorVector = pupilData.(p.Results.fitLabel).ellipses.RMSE;
                binCounts = histcounts(errorVector,histEdges);
                h=area(histEdges(1:end-1),binCounts);
                h.FaceColor = [0.5 0 0];
                hold on
                h=plot([histEdges(end) histEdges(end)],[0 sum(isnan(errorVector))],'LineWidth',2);
                set(h,'Color',[0.5 0 0]);
                plot(histEdges(end), sum(isnan(errorVector)),'*',...
                'MarkerEdgeColor',[0.5 0 0],...
                'MarkerFaceColor',[0.5 0 0]);
                % Set the 'visible' property 'off'
                ax = gca;
                ax.Visible = 'off';
                outerpos = ax.OuterPosition;
                ti = ax.TightInset;
                left = outerpos(1) + ti(1);
                bottom = outerpos(2) + ti(2);
                ax_width = outerpos(3) - ti(1) - ti(3);
                ax_height = outerpos(4) - ti(2) - ti(4);
                ax.Position = [left bottom ax_width ax_height];
                [~, pupilFileName] = fileparts(pupilFullFileName);
                text(0,0.8,pupilFileName,'Interpreter', 'none','HorizontalAlignment','left','VerticalAlignment','middle','Units','normalized');
                acqStats=sprintf('%d frames, %.1f%% nans',length(errorVector),100*sum(isnan(errorVector))/length(errorVector));
                text(0.95,0.8,acqStats,'Interpreter', 'none','HorizontalAlignment','right','VerticalAlignment','middle','Units','normalized');
                if isempty(plotPos)
                    plotPos = get(sp, 'Position');
                else
                    tmpPos = get(sp, 'Position');
                    tmpPos([1 3 4])=plotPos([1 3 4]);
                    set(sp, 'Position',tmpPos);                    
                end
            end
        end
        sp = subplot(p.Results.numPlotRows+1,1,p.Results.numPlotRows+1,'align');
        plot([histEdges(1) histEdges(end-1)],[0 0],'.w')
        hold on
        plot(histEdges(end),0,'*k')
        xlabel(['RMSE - ' p.Results.fitLabel]);
        ax = gca;
        outerpos = ax.OuterPosition;
        ti = ax.TightInset;
        left = outerpos(1) + ti(1);
        bottom = outerpos(2) + ti(2);
        ax_width = outerpos(3) - ti(1) - ti(3);
        ax_height = outerpos(4) - ti(2) - ti(4);
        ax.Position = [left bottom ax_width ax_height];
        box off
        set(gca,'YTickLabel',[],'ytick',[],'color','none');
                            tmpPos = get(sp, 'Position');
                    tmpPos([1 3 4])=plotPos([1 3 4]);
                    set(sp, 'Position',tmpPos);                    

        suptitleTT(uniqueDirNames{ii})
        %% Save the plot
        plotFileName =fullfile(plotSaveDir,[nameTags{ii} '_' p.Results.fitLabel '.pdf']);
        saveas(figHandle,plotFileName)
        close(figHandle)
    end    
end


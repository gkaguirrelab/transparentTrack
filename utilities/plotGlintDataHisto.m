function plotGlintDataHisto( dataRootDir, plotSaveDir, varargin )
% Saves plots of the quality of ellipse fitting to the glint perimeter
%
% Syntax:
%  plotGlintDataHisto( dataRootDir, plotSaveDir )
%
% Description:
%   Given a directory, the routine will perform a recursive search through
%   sub-directories to find all files with the suffix "_glint.mat".
%   Multiple glint data files within a sub-directory will be treated as
%   different acquisitions from a single session. A plot is generated for
%   each session that presents measures of the quality of the initial
%   ellipse fit to the glint perimeter for each acquisition. Each session
%   results in a separate plot, all of which are saved in the directory
%   specified by "plotSaveDir". If the plotSaveDir variable is empty, the
%   plots will be displayed and not saved.
%
%   Included local functions are subdir, by Kelly Kearney, and suptitle,
%   which is a MATLAB helper function.
%
% Inputs:
%   dataRootDir           - Full path to a directory that contains glint
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
%  'numMaxAcqs'           - The number of acquisitions displayed on a plot.
%                           This should be set equal to the maximum number
%                           of acquisitions that are expected across all
%                           sessions.
%  'fitLabel'             - The field of glintData for which the RMSE
%                           ellipse fit values are to be analyzed. The
%                           default value is for the initial ellipse fit.
%
% Outputs:
%   none
%
% Examples:
%{
    % ETTBSkip -- This is an interactive example.
    plotGlintDataHisto('','glintDataQAPlots')
%}

%% input parser
p = inputParser;

% Required
p.addRequired('dataRootDir',@ischar);
p.addOptional('plotSaveDir',[],@(x)(isempty(x) || ischar(x)));

% parse
p.parse(dataRootDir,plotSaveDir,varargin{:})


% Obtain the paths to all of the glint data files within the specified
% directory, including within sub-drectories.
fileListStruct=dir(fullfile(dataRootDir,'**','*fMRI*_glint.mat'));

% If we found at least one glint data file, then proceed.
if ~isempty(fileListStruct)
    
    % Get a list of the directories that contain a glint data file. We will
    % treat each of these as a session, and the glint data files within as
    % different acquisitions from that session.
    dirListCell=extractfield(fileListStruct,'folder');
    uniqueDirNames=unique(dirListCell);
    
    % Create a list of name tags that are the unique character components
    % of the uniqueDirNames. We trim off all characters that are shared
    % across the full path of the dir names.
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
    
    % Replace file system delimeters with valid characters. That way, we
    % can use these name tags as filenames for the resulting plots.
    nameTags=cellfun(@(x) strrep(x,filesep,'_'),nameTags,'UniformOutput',false);
    
    % Intialize a variable that will match the right-hand position of all
    % sub-plots across all plots.
    plotPos = [];
    
    % Make a temp figure
    tmpFig = figure('Position', [10 10 640 480],'visible','off');
    set(gcf,'Color',[1 1 1]);
        
    % Loop through the set of sessions
    for ii = 1:length(uniqueDirNames)
        
        % Get the list of acqusition file names
        acqList = find(strcmp(dirListCell,uniqueDirNames{ii}));
        
        % Clear the montage variable
        clear framesToMontage
        
        % Loop through the acquisitions
        for jj = 1: length(acqList)
            
            % Obtain the path to this glint data file and load it
            glintFullFileName = fullfile(fileListStruct(acqList(jj)).folder,fileListStruct(acqList(jj)).name);
            dataLoad=load(glintFullFileName);
            glintData=dataLoad.glintData;
            clear dataLoad
            
            % Obtain the modification date for this glint file
            modificationDate = glintData.meta.timestamp;
            
            % Activate the temp figure
            set(0, 'CurrentFigure', tmpFig)
            
            hold off
            scatter(glintData.X,glintData.Y,0.5,'r',...
                'filled',...
                'MarkerEdgeColor','none',...
                'MarkerFaceAlpha',0.1);
            
            set(gca,'Ydir','reverse')
            
            xlim([1 640]);
            ylim([1 480]);
            axis off
            box off
            set(gca,'LooseInset',get(gca,'TightInset'));
            gFM = glintData.meta.glintFrameMask;
            
            % Add the glint frame mask
            if ~isempty(gFM)
                recDims = [gFM(4) gFM(1) 640-gFM(2)-gFM(4) 480-gFM(3)-gFM(1)];
                hold on
                rectangle('Position',recDims)
            end
            
            % Add a blue rectangle around the entire frame
            rectangle('Position',[2 2  639 479],'EdgeColor','b','LineWidth',2);
            
            % Get the frame
            thisFrame=getframe(tmpFig);
            
            % Add a text label for the acquisition
            thisFrame.cdata = insertText(thisFrame.cdata,[20 20],fileListStruct(acqList(jj)).name,'FontSize',50);

            % Add a text label for the modification date
            thisFrame.cdata = insertText(thisFrame.cdata,[20 100],modificationDate,'FontSize',40);
            
            framesToMontage(:,:,:,jj) = imresize(thisFrame.cdata,[480 640]);
            
        end % loop over acqusitiions
        
        
        % Prepare a figure. If a plotSaveDir was not specified, then
        % display the figures. Otherwise, hide them.
        if isempty(plotSaveDir)
            figHandle=figure('visible','on');
        else
            figHandle=figure('visible','off');
        end
        
        % Format the figure
        set(gcf,'Color',[1 1 1]);
        set(figHandle, 'Units','inches')        
        height = 11;
        width = 6;
        set(figHandle, 'Position',[25 5 width height],...
            'PaperSize',[width height],...
            'PaperPositionMode','auto',...
            'Color','w',...
            'Renderer','painters'...
            );
        figure(figHandle)
        montage(framesToMontage);
        
        % Add a super-title to the entire plot
        suptitleTT(uniqueDirNames{ii})
        
        % Save the plot if a plotSaveDir has been defined
        if ~isempty(plotSaveDir)
            plotFileName =fullfile(plotSaveDir,[nameTags{ii}  '.png']);
            print(figHandle,plotFileName,'-dpng','-r600')
            close(figHandle)
        end
        
    end
end

end % plotGlintDataQA

%%% LOCAL FUNCTIONS



function hout=suptitleTT(str)
%SUPTITLE puts a title above all subplots.
%
%	SUPTITLE('text') adds text to the top of the figure
%	above all subplots (a "super title"). Use this function
%	after all subplot commands.
%
%   SUPTITLE is a helper function for yeastdemo.

%   Copyright 2003-2014 The MathWorks, Inc.

% Warning: If the figure or axis units are non-default, this
% function will temporarily change the units.

% Parameters used to position the supertitle.

% Amount of the figure window devoted to subplots
plotregion = .92;

% Y position of title in normalized coordinates
titleypos  = .95;

% Fontsize for supertitle
fs = get(gcf,'defaultaxesfontsize')+4;

% Fudge factor to adjust y spacing between subplots
fudge=1;

haold = gca;
figunits = get(gcf,'units');

% Get the (approximate) difference between full height (plot + title
% + xlabel) and bounding rectangle.

if ~strcmp(figunits,'pixels')
    set(gcf,'units','pixels');
    pos = get(gcf,'position');
    set(gcf,'units',figunits);
else
    pos = get(gcf,'position');
end
ff = (fs-4)*1.27*5/pos(4)*fudge;

% The 5 here reflects about 3 characters of height below
% an axis and 2 above. 1.27 is pixels per point.

% Determine the bounding rectangle for all the plots

h = findobj(gcf,'Type','axes');

oldUnits = get(h, {'Units'});
if ~all(strcmp(oldUnits, 'normalized'))
    % This code is based on normalized units, so we need to temporarily
    % change the axes to normalized units.
    set(h, 'Units', 'normalized');
    cleanup = onCleanup(@()resetUnits(h, oldUnits));
end

max_y=0;
min_y=1;
oldtitle = [];
numAxes = length(h);
thePositions = zeros(numAxes,4);
for i=1:numAxes
    pos=get(h(i),'pos');
    thePositions(i,:) = pos;
    if ~strcmp(get(h(i),'Tag'),'suptitle')
        if pos(2) < min_y
            min_y=pos(2)-ff/5*3;
        end
        if pos(4)+pos(2) > max_y
            max_y=pos(4)+pos(2)+ff/5*2;
        end
    else
        oldtitle = h(i);
    end
end

if max_y > plotregion
    scale = (plotregion-min_y)/(max_y-min_y);
    for i=1:numAxes
        pos = thePositions(i,:);
        pos(2) = (pos(2)-min_y)*scale+min_y;
        pos(4) = pos(4)*scale-(1-scale)*ff/5*3;
        set(h(i),'position',pos);
    end
end

np = get(gcf,'nextplot');
set(gcf,'nextplot','add');
if ~isempty(oldtitle)
    delete(oldtitle);
end
axes('pos',[0 1 1 1],'visible','off','Tag','suptitle');
ht=text(.5,titleypos-1,str);set(ht,'horizontalalignment','center','fontsize',fs,'Interpreter', 'none');
%set(gcf,'nextplot',np);
%axes(haold);
if nargout
    hout=ht;
end
end % suptitleTT


function resetUnits(h, oldUnits)
% Reset units on axes object. Note that one of these objects could have
% been an old supertitle that has since been deleted.
valid = isgraphics(h);
set(h(valid), {'Units'}, oldUnits(valid));
end


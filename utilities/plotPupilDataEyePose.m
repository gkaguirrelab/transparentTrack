function plotPupilDataQA( dataRootDir, plotSaveDir, varargin )
% Saves plots of the quality of ellipse fitting to the pupil perimeter
%
% Syntax:
%  plotPupilDataQA( dataRootDir, plotSaveDir )
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
%   Included local functions are subdir, by Kelly Kearney, and suptitle,
%   which is a MATLAB helper function.
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
%  'reportType'           - Char vector. Determines the type of plots that
%                           are created. Options include: {'rmseHisto',
%                           'eyePose'}
%
% Outputs:
%   none
%

%% input parser
p = inputParser;

% Required
p.addRequired('dataRootDir',@ischar);
p.addOptional('plotSaveDir',[],@(x)(isempty(x) || ischar(x)));

% Optional
p.addParameter('histMax',4,@isnumeric);
p.addParameter('numPlotRows',20,@isnumeric);
p.addParameter('fitLabel','initial',@ischar);
p.addParameter('reportType','rmseHisto',@ischar);

% parse
p.parse(dataRootDir,plotSaveDir,varargin{:})


% Define the edges of the bins to be used for the histogram. The histogram
% is generated for RMSE values between 0 and histMax, with 50 bins.
histEdges = 0:p.Results.histMax/50:p.Results.histMax;

% Obtain the paths to all of the pupil data files within the specified
% directory, including within sub-drectories.
fileListStruct=subdir(fullfile(dataRootDir,'*_pupil.mat'));

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
    
    % Loop through the set of sessions
    for ii = 1:length(uniqueDirNames)
        
        % Get the list of acqusition file names
        acqList = find(strcmp(fileListCell(2,:),uniqueDirNames{ii}));
        
        % Prepare a figure. If a plotSaveDir was not specified, then
        % display the figures. Otherwise, hide them.
        if isempty(plotSaveDir)
            figHandle=figure('visible','on');
        else
            figHandle=figure('visible','off');
        end
        
        % Format the figure
        set(gcf,'PaperOrientation','landscape');
        set(figHandle, 'Units','inches')
        height = 6;
        width = 11;
        set(figHandle, 'Position',[25 5 width height],...
            'PaperSize',[width height],...
            'PaperPositionMode','auto',...
            'Color','w',...
            'Renderer','painters'...
            );
        
        % Loop through the acquisitions
        for jj = 1: length(acqList)
            
            % Obtain the path to this pupil data file and load it
            pupilFullFileName = fullfile(dataRootDir,fileListCell{1,acqList(jj)});
            dataLoad=load(pupilFullFileName);
            pupilData=dataLoad.pupilData;
            clear dataLoad
            
            % Obtain the modification date for this pupil file
            if isunix
                sysCommand = ['stat -t %x ' pupilFullFileName];
                [~,modificationDateString] = system(sysCommand);
                % This returns a cell array of 4 date stamps
                modificationDateString = extractBetween(modificationDateString,'"','"');
                % The second date stamp is the modification date
                modificationDate = modificationDateString{2};
            else
                modificationDate=[];
            end
            
            % Grab just the filename for this pupil data, omitting the
            % path. We will use this to label the plot
            [pupilFilePath, pupilFileName] = fileparts(pupilFullFileName);

            % Load the associated timebase file
            fileNameStem = strsplit(pupilFileName,'_pupil');
            timebaseFileName = fullfile(pupilFilePath,[fileNameStem{1},'_timebase.mat']);
            dataLoad=load(timebaseFileName);
            timebase=dataLoad.timebase;
            clear dataLoad            

            % Check to make sure that the file has the requested field for
            % the fit result
            if isfield(pupilData,p.Results.fitLabel)
                
                % Define the subplot for this acqusition
                sp=subplot(p.Results.numPlotRows+1,1,jj,'align');
                
                % Calculate the RMSE histogram
                errorVector = pupilData.(p.Results.fitLabel).ellipses.RMSE;
                binCounts = histcounts(errorVector,histEdges);
                
                % Plot the RMSE histogram as a filled area plot
                h=area(histEdges(1:end-1),binCounts);
                h.FaceColor = [0.5 0 0];
                hold on
                
                % Add a line and marker to indicate the number of nans
                h=plot([histEdges(end) histEdges(end)],[0 sum(isnan(errorVector))],'LineWidth',2);
                set(h,'Color',[0.5 0 0]);
                plot(histEdges(end), sum(isnan(errorVector)),'*',...
                    'MarkerEdgeColor',[0.5 0 0],...
                    'MarkerFaceColor',[0.5 0 0]);
                
                % Remove the chart junk
                ax = gca;
                ax.Visible = 'off';
                
                % Force the plots to have the same dimensions
                outerpos = ax.OuterPosition;
                ti = ax.TightInset;
                left = outerpos(1) + ti(1);
                bottom = outerpos(2) + ti(2);
                ax_width = outerpos(3) - ti(1) - ti(3);
                ax_height = outerpos(4) - ti(2) - ti(4);
                ax.Position = [left bottom ax_width ax_height];
                
                % Add text to label the
                text(0,0.8,[pupilFileName ' (' modificationDate ')'],'Interpreter', 'none','HorizontalAlignment','left','VerticalAlignment','middle','Units','normalized');
                acqStats=sprintf('%d frames, %.1f%% nans, %.1f%% > %.1f%,',length(errorVector),100*sum(isnan(errorVector))/length(errorVector),100*sum(errorVector>p.Results.histMax)/length(errorVector),p.Results.histMax);
                text(0.95,0.8,acqStats,'Interpreter', 'none','HorizontalAlignment','right','VerticalAlignment','middle','Units','normalized');
                
                % If this is our very first sub-plot on the very first
                % plot, then plotPos will be empty and we therefore define
                % it. Otherwise, adjust the right-hand edge of the sub-plot
                % so that this matches across all plots.
                if isempty(plotPos)
                    plotPos = get(sp, 'Position');
                else
                    tmpPos = get(sp, 'Position');
                    tmpPos([1 3 4])=plotPos([1 3 4]);
                    set(sp, 'Position',tmpPos);
                end
            end
        end % loop over acqusitiions
        
        % Add a final sub-plot row that gives the x-axis for all plots
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
        
        % Add a super-title to the entire plot
        suptitleTT(uniqueDirNames{ii})
        
        % Save the plot if a plotSaveDir has been defined
        if ~isempty(plotSaveDir)
            plotFileName =fullfile(plotSaveDir,[nameTags{ii} '_' p.Results.fitLabel '.pdf']);
            saveas(figHandle,plotFileName)
            close(figHandle)
        end
        
    end
end

end % plotPupilDataQA

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


function varargout = subdir(varargin)
%SUBDIR Performs a recursive file search
%
% subdir
% subdir(name)
% files = subdir(...)
%
% This function performs a recursive file search.  The input and output
% format is identical to the dir function.
%
% Input variables:
%
%   name:   pathname or filename for search, can be absolute or relative
%           and wildcards (*) are allowed.  If ommitted, the files in the
%           current working directory and its child folders are returned
%
% Output variables:
%
%   files:  m x 1 structure with the following fields:
%           name:   full filename
%           date:   modification date timestamp
%           bytes:  number of bytes allocated to the file
%           isdir:  1 if name is a directory; 0 if no
%
% Example:
%
%   >> a = subdir(fullfile(matlabroot, 'toolbox', 'matlab', '*.mat'))
%
%   a =
%
%   67x1 struct array with fields:
%       name
%       date
%       bytes
%       isdir
%
%   >> a(2)
%
%   ans =
%
%        name: '/Applications/MATLAB73/toolbox/matlab/audiovideo/chirp.mat'
%        date: '14-Mar-2004 07:31:48'
%       bytes: 25276
%       isdir: 0
%
% See also:
%
%   dir

% Copyright 2006 Kelly Kearney


%---------------------------
% Get folder and filter
%---------------------------

narginchk(0,1);
nargoutchk(0,1);

if nargin == 0
    folder = pwd;
    filter = '*';
else
    [folder, name, ext] = fileparts(varargin{1});
    if isempty(folder)
        folder = pwd;
    end
    if isempty(ext)
        if isdir(fullfile(folder, name))
            folder = fullfile(folder, name);
            filter = '*';
        else
            filter = [name ext];
        end
    else
        filter = [name ext];
    end
    if ~isdir(folder)
        error('Folder (%s) not found', folder);
    end
end

%---------------------------
% Search all folders
%---------------------------

pathstr = genpath_local(folder);
pathfolders = regexp(pathstr, pathsep, 'split');  % Same as strsplit without the error checking
pathfolders = pathfolders(~cellfun('isempty', pathfolders));  % Remove any empty cells

Files = [];
pathandfilt = fullfile(pathfolders, filter);
for ifolder = 1:length(pathandfilt)
    NewFiles = dir(pathandfilt{ifolder});
    if ~isempty(NewFiles)
        fullnames = cellfun(@(a) fullfile(pathfolders{ifolder}, a), {NewFiles.name}, 'UniformOutput', false);
        [NewFiles.name] = deal(fullnames{:});
        Files = [Files; NewFiles];
    end
end

%---------------------------
% Prune . and ..
%---------------------------

if ~isempty(Files)
    [~, ~, tail] = cellfun(@fileparts, {Files(:).name}, 'UniformOutput', false);
    dottest = cellfun(@(x) isempty(regexp(x, '\.+(\w+$)', 'once')), tail);
    Files(dottest & [Files(:).isdir]) = [];
end

%---------------------------
% Output
%---------------------------

if nargout == 0
    if ~isempty(Files)
        fprintf('\n');
        fprintf('%s\n', Files.name);
        fprintf('\n');
    end
elseif nargout == 1
    varargout{1} = Files;
end

end % subdir

function [p] = genpath_local(d)
% Modified genpath that doesn't ignore:
%     - Folders named 'private'
%     - MATLAB class folders (folder name starts with '@')
%     - MATLAB package folders (folder name starts with '+')

files = dir(d);
if isempty(files)
    return
end
p = '';  % Initialize output

% Add d to the path even if it is empty.
p = [p d pathsep];

% Set logical vector for subdirectory entries in d
isdir = logical(cat(1,files.isdir));
dirs = files(isdir);  % Select only directory entries from the current listing

for i=1:length(dirs)
    dirname = dirs(i).name;
    if    ~strcmp( dirname,'.') && ~strcmp( dirname,'..')
        p = [p genpath(fullfile(d,dirname))];  % Recursive calling of this function.
    end
end

end % genpath_local

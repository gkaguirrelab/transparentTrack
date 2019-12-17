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
%     - The gray line is the eyePose fit value for the specified
%       pupilData field (either sceneConstrained or radiusSmoothed)
%     - The red plot points are the eyePose fit values, but only for those
%       points that meet criteria for having a "good" fit. 
%     - For a given frame to be considered to have a "good" fit, the
%       following criteria must be met:
%         ? the RMSE of the fit to the pupil perimeter points (in units of
%           pixels) must be belowe the rmseThreshold value (default = 3)
%         ? the fit for the frame must not have hit the upper or lower
%           bounds on the eyePose params
%     - At the bottom of the plot, gray dots mark frames with high RMSE,
%       and blue dots mark points where the fit hit the eyePose bounds.
%
%   Included local functions are subdir, by Kelly Kearney.
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
    % ETTBSkip -- This is an idiosyncratic example.
    dataRootDir = '/Users/aguirre/Dropbox (Aguirre-Brainard Lab)/TOME_processing/session1_restAndStructure';
    plotSaveDir = '/Users/aguirre/Dropbox (Aguirre-Brainard Lab)/TOME_processing/session1_restAndStructure/pupilDataQAPlots_eyePose'
    plotPupilDataEyePose(dataRootDir,plotSaveDir);
%}
%{
    plotPupilDataEyePose( '', 'pupilDataQAPlots_eyePose_sceneConstrained_FLASH_Dec2019','fieldToPlot','sceneConstrained','acquisitionStem','tfMRI_FLASH','nColumns',2)
%}
%{
    plotPupilDataEyePose( '', 'pupilDataQAPlots_eyePose_sceneConstrained_MOVIE_Dec2019','fieldToPlot','sceneConstrained','acquisitionStem','tfMRI_MOVIE','nColumns',4)
%}
%{
    plotPupilDataEyePose( '', 'pupilDataQAPlots_eyePose_sceneConstrained_RETINO_Dec2019','fieldToPlot','sceneConstrained','acquisitionStem','tfMRI_RETINO','nColumns',4)
%}

%% input parser
p = inputParser;

% Required
p.addRequired('dataRootDir',@ischar);
p.addOptional('plotSaveDir',[],@(x)(isempty(x) || ischar(x)));

% Optional
p.addParameter('rmseThreshold',3,@isscalar);
p.addParameter('eyePoseParamsToPlot',[1 2 4],@isnumeric);
p.addParameter('yRangeIncrement',[5 5 0.25],@isnumeric);
p.addParameter('xLim',[-0.5 5.6],@isnumeric);
p.addParameter('yAxisLabels',{'azimuth [deg]','elevation [deg]','radius [mm]'},@iscell);
p.addParameter('nColumns',4,@isscalar);
p.addParameter('acquisitionStem','rfMRI_REST',@ischar);
p.addParameter('fieldToPlot','radiusSmoothed',@ischar);


% parse
p.parse(dataRootDir,plotSaveDir,varargin{:})

% A constant that we will use later
msecToMin = 1/1000/60;

% Obtain the paths to all of the pupil data files within the specified
% directory, including within sub-drectories.
fileListStruct=subdir(fullfile(dataRootDir,[p.Results.acquisitionStem '*_pupil.mat']));

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
    
    % Loop through the set of sessions
    for ii = 1:length(uniqueDirNames)
        
        % Get the list of acqusition file names
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
                pupilFullFileName = fullfile(fileListCell{1,acqList(nn)});
                dataLoad=load(pupilFullFileName);
                pupilData=dataLoad.pupilData;
                clear dataLoad
                
                if ~isfield(pupilData,p.Results.fieldToPlot)
                    continue
                else
                    highRMSE = pupilData.(p.Results.fieldToPlot).ellipses.RMSE > p.Results.rmseThreshold;
                    if isfield(pupilData.(p.Results.fieldToPlot).eyePoses,'fitAtBound')
                        fitAtBound = pupilData.(p.Results.fieldToPlot).eyePoses.fitAtBound;
                        good = logical(~highRMSE .* ~fitAtBound);
                    else
                        good = logical(~highRMSE);
                    end
                    
                    lb_byParam(nn) = floor(min(pupilData.(p.Results.fieldToPlot).eyePoses.values(good,p.Results.eyePoseParamsToPlot(mm))) ./ p.Results.yRangeIncrement(mm)).*p.Results.yRangeIncrement(mm);
                    ub_byParam(nn) = ceil(max(pupilData.(p.Results.fieldToPlot).eyePoses.values(good,p.Results.eyePoseParamsToPlot(mm))) ./ p.Results.yRangeIncrement(mm)).*p.Results.yRangeIncrement(mm);
                end
            end
            lb(mm) = nanmin(lb_byParam);
            ub(mm) = nanmax(ub_byParam);
        end
        
        % Loop through the acquisitions
        for jj = 1: length(acqList)
            
            % Obtain the path to this pupil data file and load it
            pupilFullFileName = fullfile(fileListCell{1,acqList(jj)});
            dataLoad=load(pupilFullFileName);
            pupilData=dataLoad.pupilData;
            clear dataLoad
            
            % Grab just the filename for this pupil data, omitting the
            % path. We will use this to label the plot
            [pupilFilePath, pupilFileName] = fileparts(pupilFullFileName);
            fileNameStem = strsplit(pupilFileName,'_pupil');
            fileNameStem = fileNameStem{1};
            
            % Load the associated timebase file
            timebaseFileName = fullfile(pupilFilePath,[fileNameStem,'_timebase.mat']);
            dataLoad=load(timebaseFileName);
            timebase=dataLoad.timebase;
            clear dataLoad
            
            % Check that there is a smoothed radius field; otherwise
            % continue
            if ~isfield(pupilData,p.Results.fieldToPlot)
                continue
            end
            
            % Obtain the vector of good and bad time points
            highRMSE = pupilData.(p.Results.fieldToPlot).ellipses.RMSE > p.Results.rmseThreshold;
            fitAtBound = false(size(highRMSE));
            if isfield(pupilData.(p.Results.fieldToPlot).eyePoses,'fitAtBound')
                fitAtBound = pupilData.(p.Results.fieldToPlot).eyePoses.fitAtBound;
            end

            good = logical(~highRMSE .* ~fitAtBound);

            % Loop over the 3 eyePose parameters to be plotted
            for kk=1:length(p.Results.eyePoseParamsToPlot)
                
                % Define the subplot for this acqusition
                subplot(length(p.Results.eyePoseParamsToPlot),p.Results.nColumns,(kk-1)*p.Results.nColumns+jj,'align');
                
                % Plot the time-series. Make the red fit dots transparent
                plot(timebase.values*msecToMin,pupilData.(p.Results.fieldToPlot).eyePoses.values(:,p.Results.eyePoseParamsToPlot(kk)),'-','Color',[0.85 0.85 0.85],'LineWidth',0.5);
                hold on
                hLineRed = plot(timebase.values(good)/1000/60,pupilData.(p.Results.fieldToPlot).eyePoses.values(good,p.Results.eyePoseParamsToPlot(kk)),'o','MarkerSize',1);
                drawnow
                hMarkerRed = hLineRed.MarkerHandle;
                hMarkerRed.FaceColorData = uint8(255*[1; 0; 0; 0.25]);
                hMarkerRed.FaceColorType = 'truecoloralpha';
                hMarkerRed.EdgeColorData = uint8([0; 0; 0; 0]);
                
                % Add the markers for high RMSE plot points
                lowY = lb(kk) + (ub(kk)-lb(kk))/20;
                hLineGray = plot(timebase.values(highRMSE)/1000/60,repmat(lowY,size(timebase.values(highRMSE))),'o','MarkerSize',0.75);
                drawnow
                if ~isempty(hLineGray)
                    hMarkerGray = hLineGray.MarkerHandle;
                    hMarkerGray.FaceColorData = uint8(255*[0.5; 0.5; 0.5; .5]);
                    hMarkerGray.FaceColorType = 'truecoloralpha';
                    hMarkerGray.EdgeColorData = uint8([0; 0; 0; 0]);
                end
                
                % Add the markers for at bound plot points                
                if isfield(pupilData.(p.Results.fieldToPlot).eyePoses,'fitAtBound')
                    hLineBlue = plot(timebase.values(fitAtBound)/1000/60,repmat(lowY,size(timebase.values(fitAtBound))),'o','MarkerSize',0.75);
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
                
                % Remove the chart junk
                if jj ~= 1
                    set(gca,'YColor','none','TickDir','out')
                else
                    % Add a y-axis label
                    ylabel(p.Results.yAxisLabels{kk});
                end
                set(gca,'TickDir','out')
                if kk == 1
                    title({fileNameStem,pupilData.(p.Results.fieldToPlot).meta.timestamp},'Interpreter', 'none');
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
            
            % Rotate the figure by 90 degrees clockwise, because I can't get the
            % MATLAB plotting routines to output the image how I want it.
            A = imread(plotFileName);
            A = rot90(A,3);
            imwrite(A,plotFileName);

        end
        
    end % loop over sessions
end % we have at least one session
end % plotPupilDataEyePose

%%% LOCAL FUNCTIONS




function nameOut = escapeFileCharacters(nameIn)
% Sanitize file strings to be used in system commands

nameOut = strrep(nameIn,' ','\ ');
nameOut = strrep(nameOut,'(','\(');
nameOut = strrep(nameOut,')','\)');
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

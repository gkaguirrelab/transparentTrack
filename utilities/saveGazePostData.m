function saveGazePostData( dataRootDir, dataSaveDir, varargin )
% Saves plots of the quality of ellipse fitting to the pupil perimeter
%
% Syntax:
%  saveGazePostData( dataRootDir, dataSaveDir )
%
% Description:

%
%   Included local functions are subdir, by Kelly Kearney.
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
    dataSaveDir=fullfile(userpath,'projects','movieGazeTOMEAnalysis','data');
    saveGazePostData( '', dataSaveDir,'acquisitionStem','tfMRI_MOVIE')
%}
%{
    dropboxBaseDir=fullfile(getpref('eyeTrackTOMEAnalysis','dropboxBaseDir'));
    dataRootDir=fullfile(dropboxBaseDir,'TOME_processing','session1_restAndStructure');
    dataSaveDir=fullfile(dataRootDir,'pupilDataQAPlots_eyePose_July2020');
    saveGazePostData( dataRootDir, dataSaveDir,'acquisitionStem','rfMRI_REST')
%}

%% input parser
p = inputParser;

% Required
p.addRequired('dataRootDir',@ischar);
p.addOptional('dataSaveDir',[],@(x)(isempty(x) || ischar(x)));

% Optional
p.addParameter('rmseThreshold',4,@isscalar);
p.addParameter('acquisitionStem','rfMRI_REST',@ischar);


% parse
p.parse(dataRootDir,dataSaveDir,varargin{:})


% Data window
deltaT = 1000/60;
dataWindowDurMsec = 336 * 1000 - deltaT;
nSamples = round(dataWindowDurMsec/deltaT);

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
            pupilFullFileName = fullfile(fileListCell{1,acqList(jj)});
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
            
            if ~isfield(sceneGeometry.screenPosition.poseRegParams,'R')
                sceneGeomFileName
                continue
            end
            
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
            % radius smoothed
            highRMSE = pupilData.radiusSmoothed.ellipses.RMSE > p.Results.rmseThreshold;
            fitAtBound = pupilData.radiusSmoothed.eyePoses.fitAtBound;
            goodRadiusSmoothed = logical(~highRMSE .* ~fitAtBound .* ~noGlint);
            
            % Convert the eyePose to gaze position
            f = sceneGeometry.screenPosition.poseRegParams.R * [pupilData.radiusSmoothed.eyePoses.values(:,1), pupilData.radiusSmoothed.eyePoses.values(:,2)]' + sceneGeometry.screenPosition.poseRegParams.t;
            
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
                gazeData.(fileNameStem).nameTags{1} = nameTags{ii};
            else
                gazeData.(fileNameStem).vq(end+1,:,:) = vq;
                gazeData.(fileNameStem).RMSE(end+1,:) = rmse;
                gazeData.(fileNameStem).nameTags{end+1} = nameTags{ii};
            end
            
        end % loop over acquisitions       
                
    end % loop over sessions
    
    dataFileName = fullfile(dataSaveDir,'gazeData.mat');
    save(dataFileName,'gazeData');
    

end % we have at least one session

end % saveGazePostData



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

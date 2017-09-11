function convertIntrisicToWorldCoordinates(inputFileName,outputFileName,referenceFileName,varargin)

% convertIntrisicToWorldCoordinates(inputFileName,outputFileName,referenceFileName)
%
% header
% % REF for intrinsic coordinates explaination:
% https://blogs.mathworks.com/steve/2013/08/28/introduction-to-spatial-referencing/
%% Parse input
p = inputParser; p.KeepUnmatched = true;

% Required
p.addRequired('inputFileName',@ischar);
p.addRequired('outputFileName',@ischar);
p.addRequired('referenceFileName',@ischar);

% Optional params
p.addParameter('inputFileType','pupil',@ischar);
p.addParameter('referenceFileType','perimeter',@ischar); % alternative gray
p.addParameter('coordinateSystem','worldCoordinates',@ischar);

% Optional display and I/O params
p.addParameter('verbosity','none',@ischar);

% Environment parameters
p.addParameter('tbSnapshot',[],@(x)(isempty(x) | isstruct(x)));
p.addParameter('timestamp',char(datetime('now')),@ischar);
p.addParameter('hostname',char(java.lang.System.getProperty('user.name')),@ischar);
p.addParameter('username',char(java.net.InetAddress.getLocalHost.getHostName),@ischar);

p.parse(inputFileName,outputFileName,referenceFileName, varargin{:});


%% load an image from the reference file, then delete it
dataLoad=load(referenceFileName);

if strcmp(p.Results.referenceFileType,'perimeter')
    referenceFrame = squeeze(dataLoad.perimeter.data(:,:,1));
    clear dataLoad
else
    %placeholder for difference reference files
    error('No file type other than ''perimeter'' is allowed at this moment')
end

%% map the coordinates systems

xWorldLimits = [0 size(referenceFrame,2)];
yWorldLimits = [0 size(referenceFrame,1)];
R = imref2d(size(referenceFrame),xWorldLimits,yWorldLimits);


%% load the file to convert, convert and save it out

switch p.Results.inputFileType
    case 'pupil'
        % load pupil data
        
        % apply corrections to appropriate fields
        
        % repackage data and save
        
    case 'glint'
         % load pupil data
        
        % apply corrections to appropriate fields
        
        % repackage data and save
        
    otherwise
        error ('This input file type is not available at the moment')
end
        

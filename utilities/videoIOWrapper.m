function videoObj = videoIOWrapper(videoFileName, varargin)
% Wrapper for VideoReader/Writer to handle smart sync DropBox files
%
% Syntax:
%  videoObj = videoIOWrapper(videoFileName)
%
% Description:
%   Files stored on DropBox using "smart sync" appear to the file system to
%   be available. Accessing the file causes the IO operation to pause while
%   the file is downloaded and instantiated as a local copy. While many IO
%   operations tolerate this delay, the VideoReader and videoWriter
%   functions in Matlab do not. To handle this, video files are first
%   "touched" to initiate the download process, and then an attempt is made
%   to open the file.
%
%   The wrapper is used for write operations in the situation where a file
%   may already exist that is then going to be overwritten.
%
% Inputs:
%   videoFileName         - Full path to the video file to be acted upon
%
% Outputs:
%   videoObj              - Handle to the video read object
%
% Optional key/value pair:
%   maxAttempts           - The number of times to try opening the video
%   ioAction              - Char vector. Valid values are {'read','write'}
%   indexedAVI            - Logical, default false. When the ioAction is
%                           set to "wrote", if this is set to true, an
%                           indexed video is saved.
%

%% input parser
p = inputParser;

% Required
p.addRequired('videoInFileName',@ischar);

% Optional params
p.addParameter('maxAttempts',4,@isnumeric);
p.addParameter('ioAction','read',@ischar);
p.addParameter('indexedAVI',false,@ischar);

% parse
p.parse(videoFileName, varargin{:})


%% Check if the file exists
% This allows us to clear out some easy cases and identify potential errors
fileExists = exist(videoFileName, 'file') == 2;

% The file doesn't exist and we want to write it. Just open the write obj
% and then return
if ~fileExists && strcmp(p.Results.ioAction,'write')
    if p.Results.indexedAVI
        videoObj = VideoWriter(videoFileName,'Indexed AVI');
    else
        videoObj = VideoWriter(videoFileName);
    end
    return
end

% The file doesn't exist and we want to read it. Issue an error.
if ~fileExists && strcmp(p.Results.ioAction,'read')
    error('videoIOWraper:fileDoesNotExist','Requested file for video reading does not exist');
end

%% Prepare the video object
% Touch the file. If the file is in the "online only" state within a
% DropBox "smartSync" directory, this action will cause the file to be
% downloaded and made local. The system will pause during this time. The
% only effect of this step will be to update the most recent access date of
% the file. This step is only available on unix-based operating systems
if isunix
    sanitizedFileName = replace(videoFileName,{' ','(',')'},{'\ ','\(','\)'});
    % In case the touch commmand elicits an error (for example, if the file
    % is read only), pass the error text to null so that it does not appear
    % in the console.
    suppressOutputString = ' >/dev/null 2>&1';
    sysCommand = ['touch -a ' sanitizedFileName suppressOutputString];
    
    % Place the touch and videoReader command in a try-catch loop and give
    % it three tries before giving up. A pause is placed between the
    % initial touch and the attempt to open the video object. The pause
    % lengthens by a minute with subsequent try attempts. The total time
    % for the attempt to fail is thus 6 minutes plus DropBox download
    % attempt times.
    stillTrying = true; tryAttempt = 0;
    while stillTrying
        try
            system(sysCommand);
            pause(tryAttempt*60);
            switch p.Results.ioAction
                case 'read'
                    videoObj = VideoReader(videoFileName);
                case 'write'
                    if p.Results.indexedAVI
                        videoObj = VideoWriter(videoFileName,'Indexed AVI');
                    else
                        videoObj = VideoWriter(videoFileName);
                    end
            end
            stillTrying = false;
        catch
            warning('videoIOWrapper:unableToAccessVideo','Attempt %d of %d to open video failed; retrying.',tryAttempt,p.Results.maxAttempts);
            tryAttempt = tryAttempt+1;
            stillTrying = tryAttempt<(p.Results.maxAttempts+1);
        end
    end
    if ~exist('videoObj','var')
        error('videoReadWrapper:unableToAccessVideo',['Unable to read ' videoFileName]);
    end
else
    % We are on a non-unix based system (e.g., Windows). Here is where we
    % could implement similar machinery for cloud-based DropBox files. At
    % present, the routine simply creates the video object.
    switch p.Results.ioAction
        case 'read'
            videoObj = VideoReader(videoFileName);
        case 'write'
            if p.Results.indexedAVI
                videoObj = VideoWriter(videoFileName,'Indexed AVI');
            else
                videoObj = VideoWriter(videoFileName);
            end
    end
end

end

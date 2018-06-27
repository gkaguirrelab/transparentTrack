function videoInObj = videoReadWrapper(videoInFileName, varargin)
% Wrapper for VideoReader to handle accessing smart sync DropBox files
%
% Syntax:
%  videoInObj = videoReadWrapper(videoInFileName)
%
% Description:
%   A files stored on DropBox using "smart sync" appears to the file system
%   to be available. Accessing the file causes the IO operation to pause
%   while the file is downloaded and instantiated as a local copy. While
%   many IO operations tolerate this delay, the VideoReader function in
%   Matlab does not. To handle this, video files are first "touched" to
%   initiate the download process, and then an attempt is made to open the
%   file.
%
% Inputs:
%   videoInFileName       - Full path the video file to be read
%
% Outputs:
%   videoInObj            - Handle to the video read object
%
% Optional key/value pair:
%   maxAttempts           - The number of times to try opening the video
%


%% input parser
p = inputParser;

% Required
p.addRequired('videoInFileName',@ischar);

% Optional params
p.addParameter('maxAttempts',4,@isnumeric);

% parse
p.parse(videoInFileName, varargin{:})


%% Prepare the video object
% Touch the file. If the file is in the "online only" state within a
% DropBox "smartSync" directory, this action will cause the file to be
% downloaded and made local. The system will pause during this time. The
% only effect of this step will be to update the most recent access date of
% the file. This step is only available on unix-based operating systems
if isunix
    sanitizedFileName = replace(videoInFileName,{' ','(',')'},{'\ ','\(','\)'});
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
            videoInObj = VideoReader(videoInFileName);
            stillTrying = false;
        catch
            warning('videoReadWrapper:unableToReadVideo','Attempt %d of %d to open video failed; retrying.',tryAttempt,p.Results.maxAttempts);
            tryAttempt = tryAttempt+1;
            stillTrying = tryAttempt<(p.Results.maxAttempts+1);
        end
    end
    if ~exist('videoInObj','var')
        error('videoReadWrapper:unableToReadVideo',['Unable to read ' videoInFileName]);
    end
else
    % We are on a non-unix based system (e.g., Windows). Here is where we
    % could implement similar machinery for cloud-based DropBox files. At
    % present, the routine simply creates the video object.
    videoInObj = VideoReader(videoInFileName);
end

end

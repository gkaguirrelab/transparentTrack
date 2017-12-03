function resizeAndCropVideo(inputVideoName, outputVideoName, varargin)
% Resize a crop a video
%
% Description:
%   This fuction crops and resizes a video according to the specified
%   parameters. Using the default options the routine will scale and crop a
%   VGA video to match the LiveTrack standard size video. The video is also
%   converted to gray, if requested.
%
% Outputs:
%   none
%
% Inputs:
%	videoInFileName       - Full path to the video to crop/resize
%   videoOutFileName      - Full path to the output .avi file
%
% Optional key/value pairs (display and I/O):
%  'verbosity'            - Level of verbosity. [none, full]
%
% Optional key/value pairs (flow control):
%  'nFrames'              - Analyze fewer than the total number of frames
%  'startFrame'           - Which frame to start on
%
% Optional key/value pairs (analysis):
%  'resizeVideo'          - [Y X] desired output video resolution. The
%                           default values reflect the LiveTrack output.
%  'cropVideo'            - [firstX firstY lastX lastY] position of first
%                           and last
%                           pixels to include in the crop. (keep default to get livetrack
%                           format)
%  'keepOriginalSize'     - Option to skip video resizing.
%  'convertToGray'        - if set to true (default), the video will also
%                           be converted to grayscale.
%
% Outputs:
%   None
%


%% parse input and define variables
p = inputParser; p.KeepUnmatched = true;

% Required
p.addRequired('inputVideoName',@isstr);
p.addRequired('outputVideoName',@isstr);

% Optional display and I/O params
p.addParameter('verbosity', 'none', @isstr);

% Optional flow control params
p.addParameter('nFrames',Inf,@isnumeric);
p.addParameter('startFrame',1,@isnumeric);

% Optional analysis params
p.addParameter('resizeVideo',[486 720]/2, @isnumeric);
p.addParameter('cropVideo', [1 1 319 239], @isnumeric);
p.addParameter('keepOriginalSize', false, @islogical);
p.addParameter('convertToGray',true,@islogical)

% parse
p.parse(inputVideoName,outputVideoName,varargin{:})

% define variables
resizeVideo = p.Results.resizeVideo;
cropVideo = p.Results.cropVideo;
keepOriginalSize = p.Results.keepOriginalSize;

%% Prepare Video

% load video
inObj = VideoReader(inputVideoName);

% create outputVideo object
outObj = VideoWriter(outputVideoName);
outObj.FrameRate = inObj.FrameRate;
open(outObj);

% option to manually set numFrames
if p.Results.nFrames == Inf
    nFrames = floor(inObj.Duration*inObj.FrameRate);
else
    nFrames=p.Results.nFrames;
end

% alert the user
if strcmp(p.Results.verbosity,'full')
    tic
    fprintf(['Resizing and cropping video. Started ' char(datetime('now')) '\n']);
    fprintf('| 0                      50                   100%% |\n');
    fprintf('.');
end

% Resize and crop, save
for ii = p.Results.startFrame:nFrames
    % increment the progress bar
    if strcmp(p.Results.verbosity,'full') && mod(ii,round(nFrames/50))==0
        fprintf('.');
    end
    thisFrame = readFrame(inObj);
    if p.Results.convertToGray
        thisFrame = rgb2gray(thisFrame);
    end
    if keepOriginalSize == 0
        tmp = imresize(thisFrame,resizeVideo);
        thisFrame = imcrop(tmp,cropVideo);
    end
    writeVideo(outObj,thisFrame);
end

clear inObj outObj

% report completion of analysis
if strcmp(p.Results.verbosity,'full')
    fprintf('\n');
    toc
    fprintf('\n');
end

end % function
function deinterlaceVideo (videoInFileName, videoOutFileName, varargin)
% Deinterlace NTSC DV 30Hz video
%
% Syntax:
%  deinterlaceVideo (videoInFileName, videoOutFileName)
%
% Description:
%   This function deinterlaces NTSC DV 30Hz videos, saving out progressive
%   60 Hz videos, using a "bob deinterlacing" strategy.
%
%   The video is also converted to gray scale.
%
%   Four deinterlace strategies are available (bobMode):
%       'Raw' - extract 2 fields for every frame. Save progressive video.
%               Final spatial resolution is half the original resolution.
%       'Zero' - extract 2 fields for every frame. Alternate every row with
%               a row of zeros to preserve aspect ratio.
%       'Double' - extract 2 fields for every frame. Duplicate each raw to
%               preserve aspect ratio.
%       'Mean' - extract 2 fields for every frame. Add a row with the mean
%               of two consecutive rows to preserve aspect ratio.
%
%   References on bob techniques for deinterlacing and on deinterlacing in
%   general:
%       https://www.altera.com/content/dam/altera-www/global/en_US/pdfs/literature/wp/wp-01117-hd-video-deinterlacing.pdf
%       http://www.100fps.com/
%
% Inputs:
%   videoInFileName       - Full path to the video to deinterlace
%   videoOutFileName      - Full path to the output .avi file
%
% Optional key/value pairs (display and I/O):
%  'verbose'              - Logical. Default false.
%
% Optional key/value pairs (flow control):
%  'nFrames'              - Scalar. Analyze fewer than the total number of
%                           frames
%  'startFrame'           - Scalar. The frame on which to start.
%
% Optional key/value pairs (analysis):
%  'bobMode'              - String or char. The deinterlace strategy
%  'convertToGray'        - Logical. If set to true (default), the video
%                           will also be converted to grayscale.
% Outputs:
%   None
%


%% parse input and define variables
p = inputParser; p.KeepUnmatched = true;

% Required
p.addRequired('videoInFileName',@isstr);
p.addRequired('videoOutFileName',@isstr);

% Optional display and I/O params
p.addParameter('verbose',false,@islogical);

% Optional flow control params
p.addParameter('nFrames',Inf,@isscalar);
p.addParameter('startFrame',1,@isscalar);

% Optional analysis params
p.addParameter('bobMode', 'Mean', @isstr);
p.addParameter('convertToGray',true,@islogical)

% parse
p.parse(videoInFileName,videoOutFileName,varargin{:})

% define variables
bobMode = p.Results.bobMode;

% Prepare the video object
videoInObj = videoIOWrapper(videoInFileName,'ioAction','read');

% Obtain the video parameters
if p.Results.nFrames == Inf
    nFrames = floor(videoInObj.Duration*videoInObj.FrameRate);
else
    nFrames = p.Results.nFrames;
end
    
% Create the video out object
Bob = videoIOWrapper(videoOutFileName,'ioAction','write');
Bob.FrameRate = videoInObj.FrameRate * 2;
Bob.Quality = 100;
open(Bob)

% Alert the user
if p.Results.verbose
    tic
    fprintf(['Deinterlacing video. Started ' char(datetime('now')) '\n']);
    fprintf('| 0                      50                   100%% |\n');
    fprintf('.');
end

% Loop through the frames
for ii = p.Results.startFrame:nFrames
    
    % update progressbar
    if p.Results.verbose && mod(ii,round(nFrames/50))==0
        fprintf('.');
    end
    
    % get the frame
    if ~hasFrame(videoInObj)
        break
    else
        tmp = readFrame(videoInObj);
        
        % if required, convert to gray
        if p.Results.convertToGray
            thisFrame = rgb2gray(tmp);
        else
            thisFrame = tmp;
        end
        
        %get the fields
        oddFields = thisFrame(1:2:end,:);
        evenFields = thisFrame(2:2:end,:);
        
        %deinterlace
        switch bobMode
            case 'Raw'
                % shift the even lines to avoid "jumping" from frame to frame. (i.e
                % align the two fields)
                evenFields = cat(1,zeros(1,size(evenFields,2),'like',evenFields), evenFields(1:end-1,:));
                
            case 'Zero'
                % put zero rows in
                m = 1;
                k = 1;
                n = size(oddFields);
                oddFields = reshape([reshape(oddFields,m,[]);zeros(k,n(1)/m*n(2))],[],n(2));
                evenFields = reshape([reshape(evenFields,m,[]);zeros(k,n(1)/m*n(2))],[],n(2));
                evenFields = cat(1,zeros(1,size(evenFields,2),'like',evenFields), evenFields(1:end-1,:));
                
            case 'Double'
                % duplicate each row
                oddFields = repelem(oddFields, 2, 1);
                evenFields = repelem(evenFields, 2, 1);
                evenFields = cat(1,zeros(1,size(evenFields,2),'like',evenFields), evenFields(1:end-1,:));
                
            case 'Mean'
                % put means in between rows (odd fields)
                tmp = [oddFields(1,:); ((oddFields(1,:)+oddFields(2,:))/2);oddFields(2,:)];
                for jj = 2 : size(oddFields,1)-1
                    newLines = [mean([oddFields(jj,:);oddFields(jj+1,:)],1);oddFields(jj+1,:)];
                    tmp = cat(1,tmp,newLines);
                end
                oddFields = cat(1,tmp,oddFields(end,:));
                clear tmp
                clear newLines
                % put means in between rows (even fields)
                tmp = [evenFields(1,:); ((evenFields(1,:)+evenFields(2,:))./2);evenFields(2,:)];
                for jj = 2 : size(evenFields,1)-1
                    newLines = [mean([evenFields(jj,:);evenFields(jj+1,:)],1);evenFields(jj+1,:)];
                    tmp = cat(1,tmp,newLines);
                end
                evenFields = cat(1,evenFields(1,:),tmp);
                clear tmp
                clear newLines
            otherwise
                error('Unknown bobMode. Type help deinterlaceVideo for available deinterlacing methods.')
        end
        
        % write the fields as frames
        writeVideo(Bob,oddFields);
        writeVideo(Bob,evenFields);
    end
end

% report completion of analysis
if p.Results.verbose
    fprintf('\n');
    toc
    fprintf('\n');
end

% close the output video object
clear Bob

% close the input video object
clear videoInObj

end % function

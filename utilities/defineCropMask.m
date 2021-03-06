function maskBounds = defineCropMask(grayVideoName, varargin)
% Displays a frame of a video and invites the user to snap a mask box
%
% Syntax:
%  maskBounds = defineCropMask(grayVideoName)
%
% Description:
%
%
% Inputs:
%	grayVideoName         - Full path to the video in which to track the
%
% Optional key/value pairs (flow control)
%  'startFrame'           - First frame from which to start the analysis.
%

%% parse input and define variables

p = inputParser; p.KeepUnmatched = true;

% Required
p.addRequired('grayVideoName',@isstr);

% Optional flow control params
p.addParameter('startFrame',1,@isnumeric);

% parse
p.parse(grayVideoName, varargin{:})


% Prepare the video object
videoInObj = videoIOWrapper(grayVideoName,'ioAction','read');

% Grab the image dimensions
width = videoInObj.Width;
height = videoInObj.Height;

% read the video desired frame into memory, adjust gamma
thisFrame = read(videoInObj,p.Results.startFrame);

% close the video object
clear videoInObj

% create a figure for display
figureHandle=figure();

% Show the image
imageHandle = imshow(thisFrame,'Border', 'tight', 'InitialMagnification', 200);

roi = images.roi.Rectangle(gca,'Position',[round(width/4), round(height/4),round(width/2),round(height/2)]);

% Enter a while loop
notDoneFlag = true;

% Wait until the user is done (presses return or hits 'q')
while notDoneFlag
    
    % Wait for operator input
    waitforbuttonpress
    keyChoiceValue = double(get(gcf,'CurrentCharacter'));
    if ~isempty(keyChoiceValue)
        switch keyChoiceValue
            case {13,27}
                notDoneFlag = false;
        end
    end
end

% The roi position is given as Xpos, Ypos, width, height. [Xpos, Ypos] of
% [0, 0] is the upper left of the displayed image.
x = round(roi.Position);

% Convert x to the system for mask bounds, which is:
%   [nRowsTop nColumnsRight nRowsBottom nColumnsLeft]
maskBounds = [x(2) width-(x(1)+x(3)) height-(x(2)+x(4)) x(1)];
fprintf('\n');
fprintf([grayVideoName ': mask box = [%0.0f %0.0f %0.0f %0.0f]\n'],maskBounds(1),maskBounds(2),maskBounds(3),maskBounds(4));

% Clean up
close(figureHandle)


end % main function

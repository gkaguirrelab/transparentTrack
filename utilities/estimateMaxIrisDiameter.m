function maxIrisDiamPixels = estimateMaxIrisDiameter(varargin)
% Click on the an image of the eye to measure maximum iris diameter 
%
% Syntax:
%  maxIrisDiamPixels = estimateMaxIrisDiameter(grayVideoName)
%
% Description:

%
% Inputs:
%   grayVideoName         - String vector. The full path to a video file
%                           for which parameters are to be derived. If an
%                           empty string is passed (''), the operator can
%                           choose the relevant video file (mp4 or mov) via
%                           a dialog box.
%
% Outputs:
%   maxIrisDiamPixels     - Scalar
%

%% Input parser
p = inputParser;

% Required
p.addOptional('grayVideoName', [], @(x)(isempty(x) || ischar(x)));

% parse
p.parse(varargin{:})

if isempty(p.Results.grayVideoName)
    [fileName, path] = uigetfile({'*.mp4;*.mov'});
    grayVideoName = [path, fileName];
else
     grayVideoName = p.Results.grayVideoName;
end


%% Open the video 
close all
videoInObj = videoIOWrapper(grayVideoName,'ioAction','read');

%% Select and display a frame
irisDiameterFrame = GetWithDefault('>> Enter the frame number in which the visible iris diameter is largest.', [1]);
videoInObj.CurrentTime = (irisDiameterFrame - 1)/(videoInObj.FrameRate);
thisFrameIris = readFrame(videoInObj);
thisFrameIris = rgb2gray(thisFrameIris);
thisFrameIris = squeeze(thisFrameIris);
figure;
imshow(thisFrameIris, 'Border', 'tight');
hold on

%% Click on the image
fprintf('Define the maximum visible iris diameter in the figure.\n')
string = sprintf('Define the iris diameter by clicking twice on the outer boundary of the iris.');
hText = text(1,10,string, 'FontSize', 16, 'BackgroundColor', 'white');
[x1,y1] = ginput(1);
iHandle = plot(x1 ,y1, '+', 'Color', 'red');
[x2,y2] = ginput(1);
iHandle = plot(x2 ,y2, '+', 'Color', 'red');

maxIrisDiamPixels = round(sqrt((x1-x2)^2+(y1-y2)^2));

fprintf(['maxIrisDiamPixels ' num2str(maxIrisDiamPixels) ', ' grayVideoName '\n']);

end % MAIN


%%% LOCAL FUNCTIONS

function inputVal = GetWithDefault(prompt,defaultVal)
% inputVal = GetWithDefault(prompt,defaultVal)
%
% Prompt for a number or string, with a default returned if user
% hits return.
%
% 4/3/10  dhb  Wrote it.

if (ischar(defaultVal))
    inputVal = input(sprintf([prompt ' [%s]: '],defaultVal),'s');
else
    inputVal = input(sprintf([prompt ' [%g]: '],defaultVal));
end
if (isempty(inputVal))
    inputVal = defaultVal;
end

end % GetWithDefault

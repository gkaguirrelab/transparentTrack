function correctPupilPerimeterVideo(perimeterVideo,controlFileName,glintFile, correctedPerimeterVideo, varargin)

% correctPupilPerimeterVideo applies the instructions from the control file
% on the pupil perimeter video. A new corrected perimeter video will be
% saved out in the specified file.
% 
% Each frame of the original perimeter video is loaded and elaborated
% according to the control file instructions (if any) for the given frame.
% 
% Currently available instructions include:
% - there is a blink >> save out a black frame
% - there is a manually set ellipse >> draw an ellipse according to those
%       params and save the frame
% - there are cutting instruction >> cut the perimeter according to current 
%       glint position and cut instruction using the function cutPupil.m
% 
% Note that each line of the control file is set of instructions for one
% specifical video frame, identified by the FrameNumber. If there is no
% instruction for a given frame, the frame will be saved as is. The control
% file may contain multiple instruction lines referred to the same frame
% (e.g. first do a vertical cut, then do an horizontal cut); in this case,
% the routine will process the instruction on the frame one after the other
% in the order they are presented.
% 
% Once the instructions have been applied the video will be saved out.
% 
% 
% Output
% ======
%   correctedPerimeterVideo (saved where specified)
% 
% Input params
% ============
%    perimeterVideo :path to the avi perimeter video.
%    controlFileName : path to the control file (with/without extestion)
%    glintFile : path to the glint file
%    correctedPerimeterVideo : savepath for the pupil file that will be created
%    
% Usage example
% =============
%
%  correctPupilPerimeterVideo(perimeterVideo,controlFileName,glintFile, correctedPerimeterVideo)

%% parse input and define variables

p = inputParser;
% add keep unmatched option for bayes params
p.KeepUnmatched = true;
% required input
p.addRequired('perimeterVideo',@isstr);
p.addRequired('controlFileName',@isstr);
p.addRequired('glintFile',@isstr);
p.addRequired('correctedPerimeterVideo',@isstr);

%parse
p.parse(perimeterVideo, controlFileName, glintFile,correctedPerimeterVideo, varargin{:})

%% load and prepare data

% load pupilPerimeter
inObj = VideoReader(perimeterVideo);

% get number of frames
numFrames = floor(inObj.Duration*inObj.FrameRate);

% set output video object
outObj = VideoWriter(correctedPerimeterVideo);
outObj.FrameRate = inObj.FrameRate;
open(outObj);

% check controlFileName format
[~,~,ext] = fileparts(controlFileName);
if isempty(ext)
    controlFileName = [controlFileName '.mat'];
    load(controlFileName)
elseif strcmp(ext,'.mat')
    load(controlFileName)
elseif strcmp(ext,'.csv')
    controlTable  = readtable(controlFileName);
else
    error (' Only mat or csv estension can be used')
end

% load control file
controlStruct = table2struct(controlTable);

% load glint file
load(glintFile)


%% do frame by frame tracking following control file instructions
% initialize progress bar
progBar = ProgressBar(numFrames,'Correcting pupil perimeter video...');

% loop through video frames
for ff = 1:numFrames
    
    % read the frame
    thisFrame = readFrame(inObj);

    % convert to gray
    thisFrame = rgb2gray (thisFrame);
    % chech if any instruction about this frame exists
    instructionLines = find ([controlStruct.Frame] == ff);
    
    if isempty (instructionLines) % check control file
        
        % just write the frame to the correctedPerimeterVideo
        writeVideo(outObj,thisFrame);
        
         % increment progress bar
        if ~mod(ff,10);progBar(ff);end;
    
    else
        % convert frame to image
        
        img = thisFrame;
        
        % read and apply all instruction lines
        il = 1; % while counter
        while il <= length(instructionLines)
            
            % check if blink
            if controlStruct(instructionLines(il)).isBlink
                % create a black frame
                img = zeros(size(img));
                % write frame
                thisFrame = im2uint8(img);
                writeVideo(outObj,thisFrame);
                % exit while loop with flag
                doneWithThis = 1;
                break
            end %check if blink
            
            % check if force ellipse
            if ~isnan (controlStruct(instructionLines(il)).ForceEllipse_1)
                % start from back frame
                img = zeros(size(img));
                % get the ellipse params
                N = 100;
                cx = controlStruct(instructionLines(il)).ForceEllipse_1;
                cy = controlStruct(instructionLines(il)).ForceEllipse_2;
                a = controlStruct(instructionLines(il)).ForceEllipse_3;
                b = controlStruct(instructionLines(il)).ForceEllipse_4;
                phi = controlStruct(instructionLines(il)).ForceEllipse_5;
                % find ellipse points
                [Xe,Ye] = ellipse(N, cx, cy, a, b, phi);
                Xe = round(Xe);
                Ye = round(Ye);
                % draw ellipse in frame
                img(sub2ind(size(img),Ye(:),Xe(:))) = 1;
                % write frame
                thisFrame = im2uint8(img);
                writeVideo(outObj,thisFrame);
                % exit while loop with flag
                doneWithThis = 1;
                break
            end %check if force ellipse
            
            % APPLY CUT
            
            % get glint position for this frame
            Xg = glint.X(ff);
            Yg = glint.Y(ff);
            
            % get cut instructions
            U = controlStruct(instructionLines(il)).U;
            R = controlStruct(instructionLines(il)).R;
            
            % transform frame in image
            binP = imbinarize(img);
            
            % cut
            [img] = cutPupil (binP,U,R,Xg,Yg);
            
            % increment while counter
            il = il+1;
            
        end % read all instuction lines
        
        % check while loop flag
        if doneWithThis
            doneWithThis = 0;
            % increment progress bar
            if ~mod(ff,10);progBar(ff);end;
            continue
        end
        
        % write frame
        thisFrame = im2uint8(img);
        writeVideo(outObj,thisFrame);
        
        % increment progress bar
        if ~mod(ff,10);progBar(ff);end;
    end % check control file
end % loop through frames

% clear video object
clear inObj outObj


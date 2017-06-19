function [pupil, pupilTrackingParams] = mainPupilRoutine(perimeterVideo,controlFileName,glintFile, pupilFile, varargin)

% mainPupilRoutine is the core routine for ellipse fitting on the pupil
% perimeter. 
% 
% Each frame of the perimeter video is loaded and elaborated
% according to the control file instructions (if any) for the given frame.
% 
% Currently available instructions include:
% - there is a blink (do not attempt to track this frame)
% - there is a manually set ellipse (do not track, just store the hardcoded
%       ellipse params)
% - there are cutting instruction (before attempting to fit cut the
%       perimeter according to current glint position and cut instruction 
%       using the function cutPupil.m)
% 
% Note that each line of the control file is set of instructions
% for one specifical video frame, identified by the FrameNumber. The
% control file may contain multiple instruction lines referred to the same
% frame (e.g. first do a vertical cut, then do an horizontal cut); in this
% case, the routine will process the instruction on the frame one after the
% other in the order they are presented.
% 
% Once the instructions have been applied, the routine will fit an ellipse
% to the surviving perimeter pixels using a bayesian approach. ===> WRITE
% MORE DETAILS ABOUT THIS PART.
% 
% 
% Output
% ======
%    pupil : struct containing track data with the fields 
%                 X = center position of pupil on x axis [px]
%                 Y = center position of pupil on y axis [px]
%                 majorAx = length of the fitted ellipse major axis [px]
%                 minorAx = length of the fitted ellipse minor axis [px]
%                 explicitEllipseParams = explicit params for the fitted
%                   ellipse
%    pupilTrackingParams : struct containing all params used
% 
% Input params
% ============
%    perimeterVideo :path to the avi perimeter video.
%    controlFileName : path to the control file (with/without extestion)
%    glintFile : path to the glint file
%    pupilFile : savepath for the pupil file that will be created
%    
%    ADD EXTRA INPUTS FOR BAYES FITTING ROUTINE
% 
% 
% Usage example
% =============
%
%  [pupil, pupilTrackingParams] = mainPupilRoutine(perimeterVideo,controlFileName,glintFile, pupilFile)

%% parse input and define variables

p = inputParser;
% add keep unmatched option for bayes params
p.KeepUnmatched = true;
% required input
p.addRequired('perimeterVideo',@isstr);
p.addRequired('controlFileName',@isstr);
p.addRequired('glintFile',@isstr);
p.addRequired('pupilFile',@isstr);

%parse
p.parse(perimeterVideo, controlFileName, glintFile,pupilFile, varargin{:})

%% load and prepare data

% pupilPerimeter
inObj = VideoReader(perimeterVideo);

% get number of frames
numFrames = floor(inObj.Duration*inObj.FrameRate);

% check controlFileName format
[~,~,ext] = fileparts(controlFileName);
if isempty(ext)
    controlFileName = [controlFileName '.mat'];
    load(controlFileName)
elseif strcmp(ext,'mat')
    load(controlFileName)
elseif strcmp(ext,'cvs')
    controlTable  = readtable(controlFileName);
else
    error (' Only mat or csv estension can be used')
end

% load control file
controlStruct = table2struct(controlTable);

% load glint file
load(glintFile)

%% initialize pupil struct

% center of pupil
pupil.X = nan(numFrames,1);
pupil.Y = nan(numFrames,1);

% pupil axes
pupil.majorAx = nan(numFrames,1);
pupil.minorAx = nan(numFrames,1);

% fitted ellipse params
pupil.explicitEllipseParams= nan(numFrames,5);

% flags
pupil.fittingError = nan(numFrames,1);


%% do frame by frame tracking following control file instructions

% loop through video frames
for ff = 1:numFrames
    
    % EXTRACT VALID PERIMETER PIXELS
    
    % chech if any instruction about this frame exists
    instructionLines = find ([controlStruct.Frame] == ff);
    
    if isempty (instructionLines)
        
        % just read the frame and get it ready to be tracked
        thisFrame = read(inObj,ii);
        thisFrame = rgb2gray (thisFrame);
        binP = imbinarize(thisFrame);
        
    else
        % read all pertinent instruction lines
        il = 1;
        while il <= length(instructionLines)
            
            % check if blink
            if controlStruct(il).isBlink
                continue
            end %check if blink
            
            % check if force ellipse
            if ~isnan (controlStruct(il).ForceEllipse_1)
                % don't do the fitting, just store these params
                pupil.explicitEllipseParams(ff)= ...
                    [controlStruct(il).ForceEllipse_1 controlStruct(il).ForceEllipse_2 controlStruct(il).ForceEllipse_3 controlStruct(il).ForceEllipse_4 controlStruct(il).ForceEllipse_5] ;
                pupil.X(ff) = controlStruct(il).ForceEllipse_1;
                pupil.Y(ff) = controlStruct(il).ForceEllipse_2;
                pupil.majorAx(ff) = controlStruct(il).ForceEllipse_3;
                pupil.minorAx(ff) = controlStruct(il).ForceEllipse_4;
                continue
            end %check if force ellipse
            
            % APPLY CUT
            
            % get glint position for this frame
            Xg = glint.X(ff);
            Yg = glint.Y(ff);
            
            % get cut instructions
            U = controlStruct(il).U;
            R = controlStruct(il).R;
            
            % cut
            [binP] = cutPupil (binP,U,R,Xg,Yg);
            
            % increment while counter
            il = il+1;
            
        end % read all instuction lines
        
        [Yp, Xp] = ind2sub(size(binP),find(binP));
        
    end % extract valid perimeter pixels
    
                     %%%%%%%%%%%%%%%%%%%%%%%%%%%
                     %                         %
                     %   PLUG BAYES FIT HERE   %
                     %                         %
                     %%%%%%%%%%%%%%%%%%%%%%%%%%%
                     
    % store fitting params
    
    pupil.explicitEllipseParams(ff)= bayesOutput ; % DEFINE BAYES OUTPUT PROPERLY
    pupil.X(ff) = bayesOutput; % DEFINE BAYES OUTPUT PROPERLY
    pupil.Y(ff) = bayesOutput; % DEFINE BAYES OUTPUT PROPERLY
    pupil.majorAx(ff) = bayesOutput; % DEFINE BAYES OUTPUT PROPERLY
    pupil.minorAx(ff) = bayesOutput; % DEFINE BAYES OUTPUT PROPERLY
    
end % loop through frames

% clear video object
clear inObj

%% also output a trackingParams variable

pupilTrackingParams = p.Results;
% NOTE: this can be useful to store bayes params as well

%% save tracking data in the specified file

save(pupilFile,'pupil')


                     
     
                
                
                
                
                
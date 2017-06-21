function makeControlFile(controlFileName, framesToCut, blinkFrames)

% makeControlFile: produces a csv control file
% 
% Each line of the control file is called "instruction". 
% Each instruction has this format:
% FRAME NUMBER, INSTRUCTION TYPE, INSTRUCTION PARAMS
% 
% where:
%   FRAME NUMBER : frame on which to apply the instruction.
%   INSTRUCTION TYPE : what to do on the frame.
%   INSTRUCTION PARAMS : variable number of params necessary to execute the
%       instruction.
%   
% Here's the currently available instruction types and their required params:
% 'blink' - 0 params
% 'cut' - 2 params (U,R as described in cutPupil.m)
% 'ellipse' - 5 params (Xe, Ye, a, b, phi)
% '%' - param is a text string with any kind of comment.
% 
% Input params
% ============
%    controlFileName : path and name to the control file with csv extension 
%    framesToCut : array with guesses on how to cut frames, as output
%        by guessPupilCuts.m
%     blinkFrames: array with information on the blink frames as output
%        by findBlinks.m
% 
% 
% Usage example
% =============
% 
%  controlFileName = fullfile(outputDir,'ControlFile');
%  makeControlFile(controlFileName, framesToCut, blinkFrames)

%% add 1 to the blinks array
blinks = padarray(blinkFrames,[0 1], 1,'post');

%% add empty fields for cuts to the blink array

blinks = padarray(blinks,[0 2], nan,'post');

%% merge blink array into control file

cuts = [framesToCut(:,1) zeros(size(framesToCut,1),1) framesToCut(:,3) framesToCut(:,4)];

controlMat = union (cuts, blinks,'rows');

%% add columns for ellipse params
 controlMat = [controlMat nan(length(controlMat),5)];


%% make control file in a table

controlTable = array2table(controlMat,...
    'VariableNames',{'Frame' 'isBlink' 'U' 'R' 'ForceEllipse_1' 'ForceEllipse_2' 'ForceEllipse_3' 'ForceEllipse_4' 'ForceEllipse_5'});
%% save out control file as a table

save([ controlFileName '.mat'], 'controlTable')

%% save table as csv
writetable(controlTable,[ controlFileName '.csv']);
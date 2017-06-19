function makeControlFile(controlFileName, framesToCut, blinkFrames)

% makeControlFile: produces mat and csv control file

% the mat Control file is table with the following structure
% [frameNumber isBlink U R ForceEllipse_params]
% where
%   isBlink is a flag for blink frames
%   U,R are params for the cuts as described in cutPupil.m 
%   Force ellipse: 5 parameters of the manually drawn ellipse on the frame.
%   If the values are non NaN, they will override the subsequent fitting
%   steps.
% 
% Input params
% ============
%    controlFileName : path to the control file WITHOUT FILE EXTENSION
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
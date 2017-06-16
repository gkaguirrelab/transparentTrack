function makeControlFile(controlFileName, framesToCut, blinkFrames )

% makeControlFile: produces a mat and csv control file

% the mat Control file is an array with the following structure
% [frameNumber U R ForceEllipse_params]
% where U,R are params for the cuts described in cutPupil.m 
% if U and R are both NaN, that frame contains a blink.
% 
% Force ellipse: 5 parameters of the manually drawn ellipse on the frame.
% If the values are non NaN, they will override the subsequent fitting
% steps.

%% add NANs to the blinks array
blinks = padarray(blinkFrames,[0 2], nan,'post');

%% merge blink array into control file

cuts = [framesToCut(:,1) framesToCut(:,3) framesToCut(:,4)];

controlMat = union (cuts, blinks,'rows');

%% add a column to flag manual correction
 controlMat = [controlMat nan(length(controlMat),5)];

%% save out control file

save([ controlFileName '.mat'], 'controlMat')
%% make control file in a table

controlTable = array2table(controlMat,...
    'VariableNames',{'Frame' 'U' 'R' 'ForceEllipse_1' 'ForceEllipse_2' 'ForceEllipse_3' 'ForceEllipse_4' 'ForceEllipse_5'});

%% save table as csv
writetable(controlTable,[ controlFileName '.csv']);
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

%% write instructions for blinks
fid = fopen(controlFileName,'a');

for bb = 1 : length(blinkFrames)
    instruction = [num2str(blinkFrames(bb)) ',' 'blink'];
    fprintf(fid,'%s\n',instruction);
    clear instruction
end



%% write instructions for cuts

for cc = 1 : size(framesToCut,1)
    instruction = [num2str(framesToCut(cc,1)) ',' 'cut' ',' num2str(framesToCut(cc,3)) ',' num2str(framesToCut(cc,4))];
    fprintf(fid,'%s\n',instruction);
    clear instruction
end

fclose(fid);

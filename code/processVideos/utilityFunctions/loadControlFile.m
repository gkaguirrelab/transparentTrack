function instructions = loadControlFile(controlFileName)
% Imports a csv control file into a matlab variable.
%
% Description:
%   The control file format must be as specified in makeControlFile.m, i.e.
%   FRAME NUMBER, INSTRUCTION TYPE, INSTRUCTION PARAMS (variable number and
%   type).
% 
%   The instruction lines of the control files are converted in a struct as
%   follows:
%       instructions(ll).frame  = line ll, first column of csv file
%       instructions(ll).type   = line ll, second column of csv file
%       instructions(ll).params = line ll, remaining columns of csv file
%                                 (can be empty)
%
% Input: 
% 	controlFileName   - full path (with csv extension) to the control file
% 
% Output:
%   instructions      - struct with 3 fields
% 
% Example:
%	controlFileName = fullfile(outputDir,'ControlFile');
%   instructions = importControlFile(controlFileName);
% 


%% parse input and define variables
p = inputParser;

% required input
p.addRequired('controlFileName',@isstr);

% parse
p.parse(controlFileName)

% open the control file
controlFile = fopen(controlFileName);

% import values in a cell with textscan
instructionCell = textscan(controlFile,'%f%s%[^\n]','Delimiter',',');

% close the control file
fclose(controlFile);

% if there are instructions, organized them into a structure
if isempty(instructionCell{1})
    instructions(1).frame=[];
    instructions(1).type=[];
    instructions(1).params=[];
else
    frames = double(instructionCell{1});
    types =  cellstr(instructionCell{2});
    params = cellstr(instructionCell{3});
    for ff = 1: length(frames)
        instructions(ff).frame = frames(ff);
        instructions(ff).type =types{ff};
        instructions(ff).params = params{ff};
    end % loop over the number of frames
end % check if there are instructions

end % function

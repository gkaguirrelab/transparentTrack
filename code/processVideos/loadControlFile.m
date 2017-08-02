function instructions = loadControlFile(controlFileName)
% function loadControlFile(controlFileName)
%
% Imports a csv control file into a matlab variable.
% 
% The control file format must be as specified in makeControlFile.m, i.e.
% FRAME NUMBER, INSTRUCTION TYPE, INSTRUCTION PARAMS (variable number and
% type).
% 
% The instruction lines of the control files are converted in a struct as
% follows:
%       instructions(ll).frame = line ll, first column of csv file
%       instructions(ll).type = line ll, second column of csv file
%       instructions(ll).params = line ll, remaining columns of csv file
%               (can be empty)
%
% The routine will return a warning if a preliminary control file is being
% imported, but it will import it anyways.
% 
% 
% Output
% ======
%    instructions : struct with 3 fields
% Input
% =====
%    controlFileName : path and name to the control file with csv extension 
% 
% 
% Usage example
% =============
% 
%  controlFileName = fullfile(outputDir,'ControlFile');
%  instructions = importControlFile(controlFileName);
% 
%% parse input and define variables
p = inputParser;
% required input
p.addRequired('controlFileName',@isstr);

%parse
p.parse(controlFileName)

%% open the control file
controlFile = fopen(controlFileName);

%% import values in a cell with textscan
instructionCell = textscan(controlFile,'%f%s%[^\n]','Delimiter',',');

fclose(controlFile);
%% write values to the struct array
frames = double(instructionCell{1});
types =  cellstr(instructionCell{2});
params = cellstr(instructionCell{3});
for ff = 1: length(frames)
    instructions(ff).frame = frames(ff);
    instructions(ff).type =types{ff};
    instructions(ff).params = params{ff};
end

%% check if this is a preliminary control file
if strcmp(instructions(end).params, 'end of automatic instructions')
    warning ('A preliminary control file was imported. Perimeter correction might not be accurate.')
end


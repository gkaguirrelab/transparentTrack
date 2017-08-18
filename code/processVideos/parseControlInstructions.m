function [varargout] = parseControlInstructions(instructionLine,varargin)
% parseControlInstructions(instructionLine,varargin)
%
% parseInstructionParams will output the necessary paramers to execute a
% given instruction.


%% parse input

p = inputParser;

% required input
p.addRequired('instructionLine',@isstruct);

%parse
p.parse(instructionLine)


%% switch according to type of instruction

instructionType = instructionLine.type;

switch instructionType
    case '%'
        if displayComments
        fprintf('Comment for frame %d : \n',instructionLine.frame);
        fprintf('%s \n', instructionLine.params)
        end
    case 'blink'
        % just make sure that the param field is empty
         params = str2num(instructionLine.params);
        if ~isempty(params)
            error ('Invalid number of params for instruction type "blink" (0 expected)')
        end
    case 'bad'
        % just make sure that the param field is empty
         params = str2num(instructionLine.params);
        if ~isempty(params)
            error ('Invalid number of params for instruction type "bad" (0 expected)')
        end
    case 'cut'
        % parse the cut params
        params = str2num(instructionLine.params);
        if length(params) ~= 2
            error ('Invalid number of params for instruction type "cut" (2 expected)')
        else
            varargout{1} = params(1); % U
            varargout{2} = params(2); % R
            clear params
        end
    case 'ellipse'
        params = str2num(instructionLine.params);
        if length(params) ~= 5
            error ('Invalid number of params for instruction type "ellipse" (5 expected)')
        else
            varargout{1} = params(1); %Xc
            varargout{2} = params(2); %Yc
            varargout{3} = params(3); %A
            varargout{4} = params(4); %B
            varargout{5} = params(5); %phi
            
            clear params
        end
    case 'glintPatch'
        % parse the glintPatch params
        params = str2num(instructionLine.params);
        if length(params) ~= 3
            error ('Invalid number of params for instruction type "glintPatch" (3 expected)')
        else
            varargout{1} = params(1); % glint X
            varargout{2} = params(2); % glint Y
            varargout{3} = params(2); % glint patch radius
            clear params
        end
    case 'reset'
         % just make sure that the param field is empty
         params = str2num(instructionLine.params);
        if ~isempty(params)
            error ('Invalid number of params for instruction type "reset" (0 expected)')
        end
end

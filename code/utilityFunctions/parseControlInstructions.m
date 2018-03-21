function [varargout] = parseControlInstructions(instructionLine)
% Returns the parameters for a line of the control file
%
% Syntax:
%  [varargout] = parseControlInstructions(instructionLine)
%
% Description:
%   This routine accepts a single line of a control file and then returns a
%   variable number of output parameters that contain the instructions for
%   modifying a frame of a pupil perimeter file.
%
% Inputs:
%   instructionLine       - A struct that contains a line of elements from
%                           a control file
%
% Outputs:
%   varargout             - A variable number of output variables that 
%                           contain the parameters specified in the
%                           instructionLine
%


%% parse input
p = inputParser;

% required input
p.addRequired('instructionLine',@isstruct);

% parse
p.parse(instructionLine)


%% switch according to instruction type
switch instructionLine.type
    case '%'
        % this was a comment line; nothing to do
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
            varargout{3} = params(3); %area
            varargout{4} = params(4); %eccentricity
            varargout{5} = params(5); %theta
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
            varargout{3} = params(3); % glint patch radius
            clear params
        end
    case 'reset'
        % just make sure that the param field is empty
        params = str2num(instructionLine.params);
        if ~isempty(params)
            error ('Invalid number of params for instruction type "reset" (0 expected)')
        end
end % switch

end % function


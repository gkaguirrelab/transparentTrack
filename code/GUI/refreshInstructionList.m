function handles = refreshInstructionList(handles)
% this function refreshes the instruction list in the instruction box of
% the perimeter correction GUI

handles.instructionLines = find ([handles.instructions.frame] == handles.frameNumber);
    
    if isempty (handles.instructionLines) % check control file
        ListOfImageNames = {};
    else 
        ListOfImageNames = 1:1:length(handles.instructionLines);
        set(handles.instructionList,'value',1);
    end    
set(handles.instructionList,'string',ListOfImageNames);

% also display the selected instruction
instructionLine = handles.instructionLines(str2double(get(handles.instructionList,'String')));
if ~isempty(instructionLine)
    handles = loadInstructionParams(instructionLine, handles);
end

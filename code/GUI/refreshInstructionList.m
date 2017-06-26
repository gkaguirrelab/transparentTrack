function handles = refreshInstructionList(handles)
% this function refreshes the instruction list in the instruction box of
% the perimeter correction GUI

instructionLines = find ([handles.instructions.frame] == handles.frameNumber);
    
    if isempty (instructionLines) % check control file
        ListOfImageNames = {};
    else
        ListOfImageNames = 1:1:length(instructionLines);
    end
set(handles.instrunctionList,'string',ListOfImageNames);
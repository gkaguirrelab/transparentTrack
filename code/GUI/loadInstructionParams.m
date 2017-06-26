function handles = loadInstructionParams(instructionLine, handles)

% loads the instruction params in the perimeter correction GUI

switch handles.instructions(instructionLine).type
    case 'blink'
        
        % validate the instruction
        parseInstructionParams(handles.instructions(instructionLine))
        blink = 1;
        [U, R, cx, cy, a, b, phi] = deal([]);
        
    case 'ellipse'
        % get the instruction params
        [cx, cy, a, b, phi] = parseInstructionParams(handles.instructions(instructionLine));
        blink = 0;
        [U,R] = deal([]);
      
    case 'cut'
        [U,R] = parseInstructionParams(handles.instructions(instructionLine));
        blink = 0;
        [cx, cy, a, b, phi] = deal([]);
end

% set handles to the GUI elements
set(handles.blinkBtn,'value',blink)
set(handles.Utxt,'string', num2str(U))
set(handles.Rtxt,'string', num2str(R))
set(handles.ellipse1Txt,'string', num2str(cx))
set(handles.ellipse2Txt,'string', num2str(cy))
set(handles.ellipse3Txt,'string', num2str(a))
set(handles.ellipse4Txt,'string', num2str(b))
set(handles.ellipse5Txt,'string', num2str(phi))
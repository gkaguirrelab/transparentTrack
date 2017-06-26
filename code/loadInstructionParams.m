function handles = loadInstructionParams(instructionLine, handles)

% loads the instruction params in the perimeter correction GUI

switch instructionLine.type
    case 'blink'
        
        % validate the instruction
        parseInstructionParams(instructions(instructionLines(il)))
        blink = true;
        [U, R, cx, cy, a, b, phi] = deal([]);
        
    case 'ellipse'
        % get the instruction params
        [cx, cy, a, b, phi] = parseInstructionParams(instructions(instructionLines(il)));
        blink = false;
        [U,R] = deal([]);
      
    case 'cut'
        [U,R] = parseInstructionParams(instructions(instructionLines(il)));
        blink = false;
        [cx, cy, a, b, phi] = deal([]);
end

% set handles to the GUI elements
set(handles.blinkBtn,blink)
set(handles.Utxt, num2str(U))
set(handles.Rtxt, num2str(R))
set(handles.ellipse1Txt, num2str(cx))
set(handles.ellipse2Txt, num2str(cy))
set(handles.ellipse3Txt, num2str(a))
set(handles.ellipse4Txt, num2str(b))
set(handles.ellipse5Txt, num2str(phi))
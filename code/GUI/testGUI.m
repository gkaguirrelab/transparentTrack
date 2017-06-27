function varargout = testGUI(varargin)
% TESTGUI MATLAB code for testGUI.fig
%      TESTGUI, by itself, creates a new TESTGUI or raises the existing
%      singleton*.
%
%      H = TESTGUI returns the handle to a new TESTGUI or the handle to
%      the existing singleton*.
%
%      TESTGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TESTGUI.M with the given input arguments.
%
%      TESTGUI('Property','Value',...) creates a new TESTGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before testGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to testGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help testGUI

% Last Modified by GUIDE v2.5 26-Jun-2017 16:26:05

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @testGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @testGUI_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before testGUI is made visible.
function testGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to testGUI (see VARARGIN)

% parse vargin
% if a dropbox dir and param struct are passed in as variables, then load
% the corresponding data immediately
if nargin == 5
    dropboxDir = varargin{1};
    params = varargin{2};
    handles = loadDataInGUI(dropboxDir, params, handles);
end

% Choose default command line output for testGUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes testGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = testGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Frame number textbox
function frameNumTxt_Callback(hObject, eventdata, handles)
% hObject    handle to frameNumTxt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%return contents of frameNumTxt as a double
handles.frameNumber = str2double(get(hObject,'String'));

% Update handles structure
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function frameNumTxt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to frameNumTxt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in goToFrameBtn.
function goToFrameBtn_Callback(hObject, eventdata, handles)
% hObject    handle to goToFrameBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% show frames    
axes(handles.axes3);
imshow(squeeze(handles.perimVid(:,:,handles.frameNumber)));
axes(handles.axes4);
imshow(squeeze(handles.origVid(:,:,handles.frameNumber)));

% refresh slidebar
set(handles.slider1,'Value',handles.frameNumber)

% refresh instuction list
handles = refreshInstructionList(handles);

% Update handles structure
guidata(hObject, handles);


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
frameNumber = get(hObject,'Value');
handles.frameNumber = round(frameNumber);

% show images
axes(handles.axes3);
imshow(squeeze(handles.perimVid(:,:,handles.frameNumber)));
axes(handles.axes4);
imshow(squeeze(handles.origVid(:,:,handles.frameNumber)));
set(handles.slider1,'Value',handles.frameNumber)

% update frame number
set(handles.frameNumTxt, 'string', (num2str(handles.frameNumber)));

% refresh instuction list
handles = refreshInstructionList(handles);

% Update handles structure
guidata(hObject, handles);

 
% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



% --- Executes on selection change in instructionList.
function instructionList_Callback(hObject, eventdata, handles)
% hObject    handle to instructionList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns instructionList contents as cell array
%        contents{get(hObject,'Value')} returns selected item from instructionList
% 
% get the instruction line
instructionLine = handles.instructionLines(str2double(get(hObject,'String')));
if ~isempty(instructionLine)
    % load instruction
    handles = loadInstructionParams(instructionLine, handles);
    
    % display instruction
    handles = displayInstructionOnEyeVideo(handles);
    
else
    % set all handles to be empty
    % set handles to the GUI elements
    set(handles.blinkBtn,'value',0)
    set(handles.Utxt,'string', [])
    set(handles.Rtxt,'string', [])
    set(handles.ellipse1Txt,'string', [])
    set(handles.ellipse2Txt,'string', [])
    set(handles.ellipse3Txt,'string', [])
    set(handles.ellipse4Txt,'string', [])
    set(handles.ellipse5Txt,'string', [])
end

% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function instructionList_CreateFcn(hObject, eventdata, handles)
% hObject    handle to instructionList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in blinkBtn.
function blinkBtn_Callback(hObject, eventdata, handles)
% hObject    handle to blinkBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of blinkBtn
get(hObject, 'Value')

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in previewBtn.
function previewBtn_Callback(hObject, eventdata, handles)
% hObject    handle to previewBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% display instructions
handles = displayInstructionOnEyeVideo(handles);

% Update handles structure
guidata(hObject, handles);


function Utxt_Callback(hObject, eventdata, handles)
% hObject    handle to Utxt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Utxt as text
%        str2double(get(hObject,'String')) returns contents of Utxt as a double

get(hObject,'String');

% Update handles structure
guidata(hObject, handles)

% --- Executes during object creation, after setting all properties.
function Utxt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Utxt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Rtxt_Callback(hObject, eventdata, handles)
% hObject    handle to Rtxt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Rtxt as text
%        str2double(get(hObject,'String')) returns contents of Rtxt as a double

get(hObject,'String');

% Update handles structure
guidata(hObject, handles)

% --- Executes during object creation, after setting all properties.
function Rtxt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Rtxt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in deleteInstructionBtn.
function deleteInstructionBtn_Callback(hObject, eventdata, handles)
% hObject    handle to deleteInstructionBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function fileMenu_Callback(hObject, eventdata, handles)
% hObject    handle to fileMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function loadRun_Callback(hObject, eventdata, handles)
% hObject    handle to loadRun (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function saveControlFile_Callback(hObject, eventdata, handles)
% hObject    handle to saveControlFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function ellipse1Txt_Callback(hObject, eventdata, handles)
% hObject    handle to ellipse1Txt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

get(hObject,'String');

% Update handles structure
guidata(hObject, handles)
% --- Executes during object creation, after setting all properties.
function ellipse1Txt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ellipse1Txt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ellipse2Txt_Callback(hObject, eventdata, handles)
% hObject    handle to ellipse2Txt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

get(hObject,'String');

% Update handles structure
guidata(hObject, handles)

% --- Executes during object creation, after setting all properties.
function ellipse2Txt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ellipse2Txt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ellipse5Txt_Callback(hObject, eventdata, handles)
% hObject    handle to ellipse5Txt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

get(hObject,'String');

% Update handles structure
guidata(hObject, handles)

% --- Executes during object creation, after setting all properties.
function ellipse5Txt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ellipse5Txt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ellipse4Txt_Callback(hObject, eventdata, handles)
% hObject    handle to ellipse4Txt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

get(hObject,'String')

% Update handles structure
guidata(hObject, handles)

% --- Executes during object creation, after setting all properties.
function ellipse4Txt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ellipse4Txt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ellipse3Txt_Callback(hObject, eventdata, handles)
% hObject    handle to ellipse3Txt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

get(hObject,'String')

% Update handles structure
guidata(hObject, handles)

% --- Executes during object creation, after setting all properties.
function ellipse3Txt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ellipse3Txt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in saveInstructionBtn.
function saveInstructionBtn_Callback(hObject, eventdata, handles)
% hObject    handle to saveInstructionBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes on button press in drawEllipseBtn.
function drawEllipseBtn_Callback(hObject, eventdata, handles)
% hObject    handle to drawEllipseBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

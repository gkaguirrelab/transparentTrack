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

% Last Modified by GUIDE v2.5 23-Jun-2017 11:40:10

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

% load in the videos
disp('Loading videos...')
vid1 = VideoReader('~/Desktop/eyeTrackingDEMO/TOME_processing/session2_spatialStimuli/TOME_3020/050517/EyeTracking/tfMRI_FLASH_AP_run01_60Hz.avi');
vid2 = VideoReader('~/Desktop/eyeTrackingDEMO/TOME_processing/session2_spatialStimuli/TOME_3020/050517/EyeTracking/tfMRI_FLASH_AP_run01_perimeter.avi');
handles.numFrames = 100;%floor(vid2.Duration*vid2.FrameRate);
% inizialize frame array
handles.origVid = zeros([240 320 handles.numFrames],'uint8');
handles.perimVid= zeros([240 320 handles.numFrames],'uint8');
% get the videos in memory
for ii = 1:handles.numFrames
    thisFrame = readFrame(vid1);
    tmp = rgb2gray(thisFrame);
    tmp2 = imresize(tmp,[486 720]/2);
    tmp = imcrop(tmp2,[1 1 319 239]);
    handles.origVid(:,:,ii) = tmp;
end
for ii = 1:handles.numFrames
    thisFrame = readFrame(vid2);
    tmp = rgb2gray(thisFrame);
    handles.perimVid(:,:,ii) = tmp;
end
clear vid1
clear vid2
clear tmp

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


% --- Executes on button press in playbtn.
function playbtn_Callback(hObject, eventdata, handles)
% hObject    handle to playbtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pauseBtn.
function pauseBtn_Callback(hObject, eventdata, handles)
% hObject    handle to pauseBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in stopBtn.
function stopBtn_Callback(hObject, eventdata, handles)
% hObject    handle to stopBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


function frameNumTxt_Callback(hObject, eventdata, handles)
% hObject    handle to frameNumTxt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of frameNumTxt as text
%        str2double(get(hObject,'String')) returns contents of frameNumTxt as a double

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
    axes(handles.axes3);
    imshow(squeeze(handles.perimVid(:,:,handles.frameNumber)));
    axes(handles.axes4);
    imshow(squeeze(handles.origVid(:,:,handles.frameNumber)));


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


% --- Executes on button press in loadInstructionBtn.
function loadInstructionBtn_Callback(hObject, eventdata, handles)
% hObject    handle to loadInstructionBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in blinkBtn.
function blinkBtn_Callback(hObject, eventdata, handles)
% hObject    handle to blinkBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of blinkBtn


% --- Executes on button press in previewBtn.
function previewBtn_Callback(hObject, eventdata, handles)
% hObject    handle to previewBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function Utxt_Callback(hObject, eventdata, handles)
% hObject    handle to Utxt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Utxt as text
%        str2double(get(hObject,'String')) returns contents of Utxt as a double


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



function rTxt_Callback(hObject, eventdata, handles)
% hObject    handle to rTxt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of rTxt as text
%        str2double(get(hObject,'String')) returns contents of rTxt as a double


% --- Executes during object creation, after setting all properties.
function rTxt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rTxt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
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

% Hints: get(hObject,'String') returns contents of ellipse1Txt as text
%        str2double(get(hObject,'String')) returns contents of ellipse1Txt as a double


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

% Hints: get(hObject,'String') returns contents of ellipse2Txt as text
%        str2double(get(hObject,'String')) returns contents of ellipse2Txt as a double


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

% Hints: get(hObject,'String') returns contents of ellipse5Txt as text
%        str2double(get(hObject,'String')) returns contents of ellipse5Txt as a double


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

% Hints: get(hObject,'String') returns contents of ellipse4Txt as text
%        str2double(get(hObject,'String')) returns contents of ellipse4Txt as a double


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

% Hints: get(hObject,'String') returns contents of ellipse3Txt as text
%        str2double(get(hObject,'String')) returns contents of ellipse3Txt as a double


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

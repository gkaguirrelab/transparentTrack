function handles = loadDataInGUI(dropboxDir, params, handles)
% loads the requested run in the perimeter correction GUI

%% files naming
controlFileName = fullfile(dropboxDir,params.outputDir, params.projectSubfolder, ...
    params.subjectName,params.sessionDate,params.eyeTrackingDir, ...
    [params.runName '_controlFile.csv']);
glintFileName = fullfile(dropboxDir,params.outputDir, params.projectSubfolder, ...
    params.subjectName,params.sessionDate,params.eyeTrackingDir, ...
    [params.runName '_glint.mat']);
perimeterVideoName = fullfile(dropboxDir,params.outputDir, params.projectSubfolder, ...
    params.subjectName,params.sessionDate,params.eyeTrackingDir, ...
    [params.runName '_perimeter.avi']);

inputVideo = fullfile(dropboxDir,params.outputDir, params.projectSubfolder, ...
    params.subjectName,params.sessionDate,params.eyeTrackingDir, ...
    [params.runName '_60hz.avi']);
%% files loading
disp('Loading Control File')
% import control file
handles.instructions = importControlFile(controlFileName);

% import glint file (for cuts)
disp('Loading Glint File')
handles.glintFile = load(glintFileName);

% load in the videos
disp('Loading videos...')
vid1 = VideoReader(inputVideo);
vid2 = VideoReader(perimeterVideoName);
handles.numFrames = 1000;%floor(vid2.Duration*vid2.FrameRate);
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

%% first frame properties
handles.frameNumber = 1;

% set slider properties
set(handles.slider1,'Value',1,'Min',1,'Max',handles.numFrames, 'SliderStep',[1/handles.numFrames, 10/handles.numFrames])

% set frameNumberTxt Property
set(handles.frameNumTxt, 'string', (num2str(handles.frameNumber)));

% display first frame
handles.frameNumber = 1;
axes(handles.axes3);
imshow(squeeze(handles.perimVid(:,:,handles.frameNumber)));
axes(handles.axes4);
imshow(squeeze(handles.origVid(:,:,handles.frameNumber)));

% display instruction list
handles = refreshInstructionList(handles);

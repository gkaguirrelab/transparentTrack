function handles = displayInstructionOnEyeVideo(handles)

% this functions overlays the instruction information on the eye video

% blink
isblink = get(handles.blinkBtn,'value');
U = str2double(get(handles.Utxt,'string'));
R = str2double(get(handles.Rtxt,'string'));
cx = str2double(get(handles.ellipse1Txt,'string'));
cy = str2double(get(handles.ellipse2Txt,'string'));
a = str2double(get(handles.ellipse3Txt,'string'));
b = str2double(get(handles.ellipse4Txt,'string'));
phi = str2double(get(handles.ellipse5Txt,'string'));


% no instruction
if isempty(find(~isnan([U R cx cy a b phi]))) && ~isblink
    
    % plot full perim
    axes(handles.axes3);
    imshow(squeeze(handles.perimVid(:,:,handles.frameNumber)), 'Border', 'tight')
    
    % plot original frame
    axes(handles.axes4);
    imshow(squeeze(handles.origVid(:,:,handles.frameNumber)))
end

% there's a blink
if isblink && isempty(find(~isnan([U R cx cy a b phi])))
    %perim video : make frame black
    axes(handles.axes3);
    I = squeeze(handles.origVid(:,:,handles.frameNumber));
    blinkI = zeros(size(I));
    imshow(blinkI, 'Border', 'tight')
    
    
    %orig video : add blink flag
    axes(handles.axes4);
    position = [100 100];
    text_str = 'BLINK';
    RGB = insertText(squeeze(handles.origVid(:,:,handles.frameNumber)),...
        position,text_str,'FontSize',18,'BoxColor',...
        'red','BoxOpacity',0.4,'TextColor','white');
    imshow(RGB, 'Border', 'tight')
end

% there are cutting instructions
if  ~isnan(U) && ~isnan(R) && isempty(find(~isnan([cx cy a b phi]))) && ~isblink
    % get glint position
    Xg = handles.glintFile.glintData.X(handles.frameNumber);
    Yg = handles.glintFile.glintData.Y(handles.frameNumber);
    % get perimeter frame
    binP = imbinarize(squeeze(handles.perimVid(:,:,handles.frameNumber)));
    
    % cut the frame
    [binPcut] = cutPupil (binP,U,R,Xg,Yg);
    
    % perim video : display the cut
    axes(handles.axes3);
    imshow(im2uint8(binPcut), 'Border', 'tight')
    
    % orig video : overlay the cut
    axes(handles.axes4);
    imshow(squeeze(handles.origVid(:,:,handles.frameNumber)))
    green = cat(3, zeros(size(squeeze(handles.perimVid(:,:,handles.frameNumber)))), ones(size(squeeze(handles.perimVid(:,:,handles.frameNumber)))), zeros(size(squeeze(handles.perimVid(:,:,handles.frameNumber)))));
    hold on
    h = imshow(green);
    hold off
    % h = imshow(im2uint8(binPcut), 'Border', 'tight');
    set(h, 'AlphaData', binPcut);
end

% there are ellipse params
if isempty(find(isnan([cx cy a b phi]))) && isempty(find(~isnan([U R]))) && ~isblink
    % plot ellipse
    
    % plot ellipse
    
end

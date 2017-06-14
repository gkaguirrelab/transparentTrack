function showBlinks(blinkFrames,grayI)

% this functions just shows the blink flag in form of a text box over a
% given video

% get number of frames from grayI
numFrames = size(grayI,3);

ih = figure;

% loop through gray frames
for i = 1:numFrames
    
    % Get the frame
    I = squeeze(grayI(:,:,i));
    
    % Show the frame
    imshow(I, 'Border', 'tight');
    
    if ismember(i,blinkFrames)
        
        position = [100 100];
        text_str = 'BLINK';
        
        RGB = insertText(I,position,text_str,'FontSize',18,'BoxColor',...
            'red','BoxOpacity',0.4,'TextColor','white');
        imshow(RGB, 'Border', 'tight')
    end
end

close (ih)

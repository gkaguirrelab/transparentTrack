function [pupil , iris, eyelid] = irisFit (params)

%% set defaults

% params for image resizing and cropping
if ~isfield(params,'imageSize')
    params.imageSize = [486 720]/2;
end
if ~isfield(params,'imageCrop')
    params.imageCrop = [1 1 319 239];
end


%params.irisRadius=210;




%% Load video
disp('Loading video file, may take a couple minutes...');
inObj                   = VideoReader(params.inVideo);
numFrames               = floor(inObj.Duration*inObj.FrameRate);
% option to overwrite numFrames (for quick testing)
if isfield(params,'forceNumFrames')
    numFrames = params.forceNumFrames;
end

% initialize gray image array
grayI                   = zeros([240 320 numFrames],'uint8');

% Convert to gray, resize, crop to livetrack size
for i = 1:numFrames
    thisFrame           = readFrame(inObj);
    tmp                 = rgb2gray(thisFrame);
    tmp2        = imresize(tmp,params.imageSize);
    tmp3 = imcrop(tmp2,params.imageCrop);
    grayI(:,:,i) = tmp3;
end

if isfield(params,'outVideo')
    outObj              = VideoWriter(params.outVideo);
    outObj.FrameRate    = inObj.FrameRate;
    open(outObj);
end

clear RGB inObj


%% initialize variables
pupil.X = nan(numFrames,1);
pupil.Y = nan(numFrames,1);
pupil.radiusX = nan(numFrames,1);
pupil.radiusY = nan(numFrames,1);
iris.X = nan(numFrames,1);
iris.Y = nan(numFrames,1);
iris.radius = nan(numFrames,1);
eyelid.mask = zeros([240 320 numFrames],'uint8');
%% track

for i = 1:numFrames
    % Get the frame
    Image = squeeze(grayI(:,:,i));
    
    % Turn off warnings for NARGCHK. This will catch up with us eventually...
    warnID='MATLAB:nargchk:deprecated';
    orig_state = warning;
    warning('off',warnID);
    
    % Set the params to empty by default. May need to provide a irisRadius (in
    % pixels) if the routine cannot find the iris boundary consistently.
    params=[];
    %params.irisRadius=210;
    
    % Perform the segmentation
    try
        [pupilCoords, irisCoords, eyelidMask ] = irisseg_main(Image, params);
    catch ME
    end
    if exist ('ME', 'var')
        fig = figure;
        colormap gray
        subplot(2,1,1)
        image(Image);
        axis equal
        axis([0 size(Image,2) 0 size(Image,1)]);
        % save as a frame
        F = getframe(fig);
        writeVideo(outObj,F);
        close(fig);
        clear ME
        continue
    end
    
    
    % Restore the warning state
    warning(orig_state);
    
    % mark up the image with an ellipse for the pupil
    if(size(Image,3)==3)
        markedImage=rgb2gray(Image);
    else
        markedImage=Image;
    end
    
    centerX=pupilCoords(1);
    centerY=pupilCoords(2);
    radiusX=pupilCoords(3);
    radiusY=pupilCoords(4);
    ringThick=3;
    [columnsInImage rowsInImage] = meshgrid(1:size(Image,2), 1:size(Image,1));
    innerRing = (rowsInImage - centerY).^2 ./ radiusY^2 ...
        + (columnsInImage - centerX).^2 ./ radiusX^2 <= 1;
    outerRing = (rowsInImage - centerY).^2 ./ (radiusY+ringThick)^2 ...
        + (columnsInImage - centerX).^2 ./ (radiusX+ringThick)^2 <= 1;
    ellipseRingPixels= ~(outerRing <= innerRing);
    markedImage(ellipseRingPixels)=254;
    
    % add the iris circle
    markedImage=insertShape(markedImage,'Circle',[irisCoords(1) irisCoords(2) irisCoords(3)],'Color','black','LineWidth',20);
    
    % Show the segmentation
    fig = figure;
    colormap gray
    subplot(2,1,1)
    image(markedImage);
    axis equal
    axis([0 size(Image,2) 0 size(Image,1)]);
    
    subplot(2,1,2)
    imagesc(eyelidMask);
    axis equal
    axis([0 size(Image,2) 0 size(Image,1)]);
    
    % save as a frame
    F = getframe(fig);
    writeVideo(outObj,F);
    close(fig);
    
    % save variables
    pupil.X(i) = pupilCoords(1);
    pupil.Y(i) = pupilCoords(2);
    pupil.radiusX(i) = pupilCoords(3);
    pupil.radiusY(i) = pupilCoords(4);
    iris.X(i) = irisCoords(1);
    iris.Y(i) = irisCoords(2);
    iris.radius(i) = irisCoords(3);
    try
    eyelid.mask(:,:,i) = eyelidMask;
    catch ME
    end
    if exist ('ME', 'var')
        clear ME
        continue
    end
end

% save full video
if isfield(params,'outVideo')
    close(outObj);
end

if isfield(params,'outMat')
    save(params.outMat,'pupil','iris','eyelid');
end


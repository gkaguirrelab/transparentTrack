function saveEyeModelMontage(obj,fileNameSuffix)

%% Obtain variables from the object
sceneGeometry = obj.sceneGeometry;
modelEyePose = obj.modelEyePose;
frameSet = obj.frameSet;
videoStemName = obj.videoStemName;
args = obj.args;

perimeter = args{1};
gazeTargets = args{4};
grayVideoName = [videoStemName '_gray.avi'];
montageFileName = [videoStemName '_sceneGeometry_eyeModelMontage' fileNameSuffix '.png'];



% Saves a montage with the model eye superimposed.

% Silence some errors that can arise during the forward projection
warningState = warning;
warning('off','projectModelEye:ellipseFitFailed');


% We will create a square grid of frames, and try to position the frames in
% the grid by their eye pose. First, figure out how many grid squares we
% have
nFrames = length(frameSet);
nGrid = ceil(sqrt(nFrames));

% Convert the gazeTargets array into an ordered list to store the frames
pos = gazeTargets ./ max(unique(abs(gazeTargets)));

% Reverse the sign of the elevation row, as positive eye movements
% correspond to the eye moving upwards in the image
pos(2,:) = -pos(2,:);

% Convet the positions into an ordered list
montageOrder = (pos(2,:)+nGrid-1).*nGrid+pos(1,:)-1;

% Sort the ellipse array list so that the frames appear in temporal order
[frameSet, sortOrder] = sort(frameSet);
modelEyePose = modelEyePose(sortOrder,:);
perimeter.data = perimeter.data(sortOrder);

% If we have a full supply of gazeTargets, use them to define the montage
% order. Otherwise, just go with the sortOrder
if all(~isnan(sum(gazeTargets)))
    montageOrder = montageOrder(sortOrder);
else
    montageOrder = 1:length(montageOrder);
end


% Check that the file exists
if exist(grayVideoName,'file') && ~isempty(frameSet)
    
    % Open the video object
    videoInObj = videoIOWrapper(grayVideoName,'ioAction','read');
    
    % Get the video properties
    videoSizeX = videoInObj.Width;
    videoSizeY = videoInObj.Height;
    
    % Define a variable to hold the selected frames
    framesToMontage = zeros(videoSizeY,videoSizeX,3,nGrid^2,'uint8');
    
    % Define a figure
    hFig = figure( 'Visible', 'off');
    hAxes = gca();
    
    % A counter for which frame of the video is currently available to read
    curFrame = 0;
    
    % Loop through the frames and keep the matching ones
    for ii = 1:length(frameSet)
        idx = frameSet(ii);
        % Read up to the desired frame
        while curFrame < idx
            sourceFrame = readFrame(videoInObj);
            curFrame = curFrame + 1;
        end
        imshow(sourceFrame,'Border', 'tight','Parent',hAxes);
        hold on
        axis off;
        % Add the pupil perimeter
        plot(perimeter.data{ii}.Xp ,perimeter.data{ii}.Yp, '.w', 'MarkerSize', 1);
        % Add the rendered eye model
        eyePose = modelEyePose(ii,:);
        if ~any(isnan(eyePose))
            renderEyePose(eyePose, sceneGeometry, 'newFigure', false, ...
                'modelEyeLabelNames', {'retina' 'irisPerimeter' 'pupilEllipse' 'cornea' 'glint_01' 'glint_02'}, ...
                'modelEyePlotColors', {'.w' '.b' '-g' '.y' 'xr' 'xr'}, ...
                'modelEyeAlpha', [0.25 0.25 0.25 0.25 1 1],...
                'modelEyeSymbolSizeScaler',1.5,...
                'showAzimuthPlane',true);
        end
        % Get the frame
        drawnow;
        thisFrame=getframe(hFig);
        % Add a text label for the frame number
        frameLabel = sprintf('frame: %d',idx);
        thisFrame.cdata = insertText(thisFrame.cdata,[20 20],frameLabel,'FontSize',30);
        % Store the frame. Detect if we have a bad or empty frame and then
        % skip if that is the case
        if all(size(squeeze(framesToMontage(:,:,:,ii)))==size(thisFrame.cdata))
            framesToMontage(:,:,:,montageOrder(ii)) = thisFrame.cdata;
        end
        % hold off
        hold off
    end
    
    % Close the temporary figure
    close(hFig);
    
    % Prepare the figure
    figHandle=figure('Visible','off');
    set(gcf,'PaperOrientation','landscape');
    set(figHandle, 'Units','inches')
    height = 6;
    width = 11;
    
    % The last two parameters of 'Position' define the figure size
    set(figHandle, 'Position',[25 5 width height],...
        'PaperSize',[width height],...
        'PaperPositionMode','auto',...
        'Color','w');
    
    % Turn off a warning that can occur during the montage step
    warningState = warning;
    warning('off','images:initSize:adjustingMag');
    
    % Create the montage
    montage(framesToMontage);
    
    % Restore the warning state
    warning(warningState);
    
    % Save the montage
    saveas(figHandle,montageFileName)
    
    % Close the figure
    close(figHandle)
    
    % Rotate the figure by 90 degrees clockwise, because I can't get the
    % MATLAB plotting routines to output the image how I want it.
    A = imread(montageFileName);
    A = rot90(A,3);
    imwrite(A,montageFileName);
    
    % close the video object
    clear videoInObj
    
end % There is a file to plot

% Restore the warning state
warning(warningState);


end % saveEyeModelMontage

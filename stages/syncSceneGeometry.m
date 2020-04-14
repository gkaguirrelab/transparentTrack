function syncSceneGeometry(videoStemNameIn, videoStemNameOut, varargin)
% Update sceneGeometry camera position to match an acquisition
%
% Syntax:
%  syncSceneGeometry(videoStemNameIn, videoStemNameOut, varargin)
%
% Description:
%   A sceneGeometry file specifies biometric properties of an eye, and the
%   spatial arrangement of the eye with respect to a camera. This routine
%   adjusts a sceneGeometry file created for one acquisition to best match
%   the measurements made in a different acqusition.
%
%   The adjustment is guided by fitting the glint and pupil within a
%   selection of frames that are distributed across eye positions and time.
%
% Inputs:
%   videoStemNameIn       - Full path to the set of files that provide the
%                           source sceneGeometry. The path should exclude
%                           the trailing "_gray.avi" that identifies the
%                           video.
%   videoStemNameIn       - Full path to the set of files for which a new
%                           sceneGeometry file is to be created. The path
%                           should exclude the trailing "_gray.avi" that
%                           identifies the video.
%
% Optional key/value pairs (display and I/O):
%  'verbose'              - Logical.
%  'saveDiagnosticPlot'   - Logical.
%  'doNotSyncSceneToItself' - Logical. Strange things could happen if the
%                           user requests that a sceneGeometry file be
%                           synced to itself. This circumstance is detected
%                           and the routine exits unless this flag is set
%                           to false.
%
% Optional key/value pairs (flow control)
%  'useParallel'          - If set to true, use the MATLAB parallel pool
%  'nWorkers'             - Specify the number of workers in the parallel
%                           pool. If undefined the default number will be
%                           used.
%
% Optional key/value pairs (fitting params):
%  'alignMethod'          - Cell or char vec. Controls how the routine 
%                           selects frames from the pupilFile acquisition
%                           to use as the fixed target to which the
%                           sceneGeometry is adjusted. Valid options are:
%                               {'gazePre','gazePost','shape'}

%
% Examples:
%{
%}


%% Parse vargin for options passed here
p = inputParser; p.KeepUnmatched = true;

% Required
p.addRequired('videoStemNameIn',@ischar);
p.addRequired('videoStemNameOut',@ischar);

% Optional display and I/O params
p.addParameter('verbose',true,@islogical);
p.addParameter('saveDiagnosticPlot',true,@islogical);
p.addParameter('doNotSyncSceneToItself',true,@islogical);

% Optional flow control params
p.addParameter('useParallel',true,@islogical);
p.addParameter('nWorkers',[],@(x)(isempty(x) || isnumeric(x)));

% Optional fitting params
p.addParameter('alignMethod','gazePre',@(x)(ischar(x) | iscell(x)));


%% Parse and check the parameters
p.parse(videoStemNameIn, videoStemNameOut, varargin{:});


alignMethod = p.Results.alignMethod;
if iscell(alignMethod)
    alignMethod = alignMethod{1};
end

%% Check if we are syncing a sceneGeometry to itself
if p.Results.doNotSyncSceneToItself && strcmp(videoStemNameIn,videoStemNameOut)
    if p.Results.verbose
        fprintf('Detected that the sceneGeometry source is the same as the target acquisition; returning.\n');
    end
    return
end


%% Load sceneGeometry for the videoStemNameIn
% This sceneGeometry provides the eye model and the parameters that
% transform eyePose into gaze position.

% Load the sceneGeometry variable into memory
dataLoad=load([videoStemNameIn '_sceneGeometry.mat']);
sceneGeometryIn=dataLoad.sceneGeometry;
clear dataLoad

% Get the camera offset point. This is where the center of the video image
% lands on the camera sensor array. We will need this later when we make
% some diagnostic images of the alignment of the in and out videos.
cameraOffsetPoint = [sceneGeometryIn.cameraIntrinsic.matrix(1,3), ...
    sceneGeometryIn.cameraIntrinsic.matrix(2,3)];


%% Find the fixation frame for sceneGeometryIn
% We assume that the sceneGeometryIn contains a frame that has been
% assigned a gazeTarget value of [0; 0].

% Which of the list of frames is the [0;0] fixation frame
idx = find((sceneGeometryIn.meta.estimateSceneParams.obj.gazeTargets(1,:)==0).*(sceneGeometryIn.meta.estimateSceneParams.obj.gazeTargets(2,:)==0));

% Store the eyePose for this frame
eyePoseFixationIn = sceneGeometryIn.meta.estimateSceneParams.obj.modelEyePose(idx,:);

% Store the pupilEllipse for this frame
pupilEllipseFixationIn = sceneGeometryIn.meta.estimateSceneParams.obj.modelPupilEllipse(idx,:);

% Find the shape of the pupil for this frame, expressed as theta and rho
% values (SEE: csaEllipseError)
rhoIn = 1-sqrt(1-pupilEllipseFixationIn(4)^2);
thetaIn = pupilEllipseFixationIn(5)*2;

% Get the scene and eye parameters from the sceneGeometryIn
model.eye.x0 = sceneGeometryIn.meta.estimateSceneParams.xEye;
model.scene.x0 = sceneGeometryIn.meta.estimateSceneParams.xScene;
eyeArgs = sceneGeometryIn.meta.estimateSceneParams.obj.setupArgs;
errorArgs = { ...
    'poseRegParams',sceneGeometryIn.meta.estimateSceneParams.obj.poseRegParams,...
    'vecRegParams',sceneGeometryIn.meta.estimateSceneParams.obj.vecRegParams};

% Get the cameraDepth from the sceneGeometryIn
cameraDepth = sceneGeometryIn.meta.estimateSceneParams.xScene(end);

% Load in the video image for this frame.
absIdx = sceneGeometryIn.meta.estimateSceneParams.obj.frameSet(idx);
videoFrameIn = makeMedianVideoImage([videoStemNameIn '_gray.avi'],'startFrame',absIdx,'nFrames',1);


%% Select frames to guide the search
% Identify a frame in the videoStemNameOut acquisition that best
% corresponds to the eye when it is looking at the fixation point. This
% could come from an explicit fixation period, or from matching the shape
% of the pupil to that observed during the videoStemNameIn fixation period.

switch alignMethod
    case 'gazePre'
        [frameSet, gazeTargets] = selectFrames.gazePre(videoStemNameOut);
    case 'gazePost'
        [frameSet, gazeTargets] = selectFrames.gazePost(videoStemNameOut);
    case 'shape'
        [frameSet, gazeTargets] = selectFrames.shape(videoStemNameOut, rhoIn, thetaIn);
    otherwise
        error('This is not a defined align method')
end

% Load in the video image for this frame.
videoFrameOut = makeMedianVideoImage([videoStemNameOut '_gray.avi'],'startFrame',frameSet(1),'nFrames',1);

% Add frames that are distributed across time and space.
[frameSetA, gazeTargetsA] = selectFrames.gridTime(videoStemNameOut);
[frameSetB, gazeTargetsB] = selectFrames.gridSpace(videoStemNameOut);
frameSet = [frameSet frameSetA frameSetB];
gazeTargets = [gazeTargets gazeTargetsA gazeTargetsB];


%% Perform the synchronization search
estimateSceneParams(videoStemNameOut, frameSet, gazeTargets, ...
    'searchStrategy','sceneSync','cameraDepth',cameraDepth,'model',model,...
    'eyeArgs',eyeArgs,'errorArgs',errorArgs, ...
    'verbose',p.Results.verbose,...
    'useParallel',p.Results.useParallel,...
    'nWorkers',p.Results.nWorkers);

% Load the sceneGeometryOut variable into memory
dataLoad=load([videoStemNameOut '_sceneGeometry.mat']);
sceneGeometryOut=dataLoad.sceneGeometry;
clear dataLoad

% Store the eyePose for the fixation frame of sceneGeometryOut
eyePoseFixationOut = sceneGeometryOut.meta.estimateSceneParams.obj.modelEyePose(1,:);


%% Create and save a diagnostic figure
if p.Results.saveDiagnosticPlot
            
    %% videoStemNameIn
    % Typically the gazeCal source

    displayImage = videoFrameIn;
    tmpFig = figure('visible','off');
    renderEyePose(eyePoseFixationIn, sceneGeometryIn, ...
        'newFigure', false, 'visible', false, ...
        'backgroundImage',displayImage, ...
        'showAzimuthPlane', true, ...
        'modelEyeLabelNames', {'retina' 'pupilEllipse' 'cornea' 'glint_01'}, ...
        'modelEyePlotColors', {'.w' '-g' '.y' 'xr'}, ...
        'modelEyeSymbolSizeScaler',1.5,...
        'modelEyeAlpha', [0.25 0.25 0.25 1]);
    text(20,30,videoStemNameIn, 'Color', 'g','Fontsize',16,'Interpreter','none');
    msg = ['frame ' num2str(absIdx)];
    addAnnotation(msg);
    hold on
    plot([size(displayImage,2)/2, size(displayImage,2)/2],[0 size(displayImage,2)],'-b');
    plot([0 size(displayImage,1)],[size(displayImage,1)/2, size(displayImage,1)/2],'-b');
    tmpFrame = getframe(gcf);
    imageSet(1) = {tmpFrame.cdata};
    close(tmpFig);
    
    
    %% videoStemNameOut
    % The acquisition for which we have new sceneGeometry
    
    displayImage = videoFrameOut;
    tmpFig = figure('visible','off');
    renderEyePose(eyePoseFixationOut, sceneGeometryOut, ...
        'newFigure', false, 'visible', false, ...
        'backgroundImage',displayImage, ...
        'showAzimuthPlane', true, ...
        'modelEyeLabelNames', {'retina' 'pupilEllipse' 'cornea' 'glint_01'}, ...
        'modelEyePlotColors', {'.w' '-g' '.y' 'xr'}, ...
        'modelEyeSymbolSizeScaler',1.5,...
        'modelEyeAlpha', [0.25 0.25 0.25 1]);
    text(20,30,videoStemNameOut, 'Color', 'r','Fontsize',16,'Interpreter','none');
    msg = ['frame ' num2str(frameSet(1))];
    addAnnotation(msg);
    % Add cross hairs
    hold on
    plot([size(displayImage,2)/2, size(displayImage,2)/2],[0 size(displayImage,2)],'-b');
    plot([0 size(displayImage,1)],[size(displayImage,1)/2, size(displayImage,1)/2],'-b');
    tmpFrame = getframe(gcf);
    imageSet(2) = {tmpFrame.cdata};
    close(tmpFig);
    
    
    %% Difference image
    
    % Transform the videoFrameIn to the space of videoFrameOut
    regParams = calcImageTransform(sceneGeometryIn,sceneGeometryOut,cameraOffsetPoint);
    videoFrameInAdj = updateFrame(videoFrameIn,regParams,cameraOffsetPoint);
    displayImage = videoFrameOut - double(videoFrameInAdj);
    tmpFig = figure('visible','off');
    imshow(displayImage,[], 'Border', 'tight');
    text(20,30,'Difference', 'Color', 'w','Fontsize',16,'Interpreter','none');
    tmpFrame = getframe(gcf);
    imageSet(3) = {tmpFrame.cdata};
    close(tmpFig);
    
    % Prepare the figure
    figHandle=figure('visible','off');
    set(gcf,'PaperOrientation','portrait');
    
    set(figHandle, 'Units','inches')
    height = 4;
    width = 12;
    
    % The last two parameters of 'Position' define the figure size
    set(figHandle, 'Position',[25 5 width height],...
        'PaperSize',[width height],...
        'PaperPositionMode','auto',...
        'Color','w',...
        'Renderer','painters'...
        );
    
    % Post the montage of the imageSet
    montage(imageSet,'Size', [1 3]);
    
    % Post the title
    pathParts = strsplit(videoStemNameIn,filesep);
    if length(pathParts)>4
        titleString = fullfile(pathParts{end-4:end-2});
    else
        titleString = videoStemNameIn;
    end
    title(titleString,'Interpreter','none')
    
    % Report the alignment method
    annotation('textbox', [0.15, .125, 0, 0], 'string', alignMethod,'FontWeight','bold','FitBoxToText','on','LineStyle','none','HorizontalAlignment','left','Interpreter','none')
        
    % Add a text summary below. If any delta fixation angle is geater than
    % 1 deg, print the message text in red to alert that this was a large
    % eye rotation change.
    deltaX = sceneGeometryOut.meta.estimateSceneParams.xScene - sceneGeometryIn.meta.estimateSceneParams.xScene;
    deltaPose = eyePoseFixationOut - eyePoseFixationIn;
    
    msg = sprintf('delta torsion [deg] = %2.1f',deltaX(3));
    annotation('textbox', [0.5, .175, 0, 0], 'string', msg,'FitBoxToText','on','LineStyle','none','HorizontalAlignment','center','Interpreter','none')
    msg = sprintf('delta translation [mm] [x; y; z] = [%2.3f; %2.3f; %2.3f]',deltaX(4:6));
    annotation('textbox', [0.5, .125, 0, 0], 'string', msg,'FitBoxToText','on','LineStyle','none','HorizontalAlignment','center','Interpreter','none')
    msg = sprintf('delta eye pose [azi, ele, tor, radius] = [%2.3f, %2.3f, %2.3f, %2.3f]',deltaPose);
    msgColor = 'black';
    if any(abs(deltaPose(1:3)) > 1.0)
        msgColor = 'red';
    end
    annotation('textbox', [0.5, .075, 0, 0], 'string', msg,'Color',msgColor,'FitBoxToText','on','LineStyle','none','HorizontalAlignment','center','Interpreter','none')
    
    % Save and close the figure
    print(figHandle,[videoStemNameOut '_sceneGeometry_sceneSync_QA.pdf'],'-dpdf','-bestfit');
    close(figHandle);
    
end



end % Main function





%% LOCAL FUNCTIONS


function regParams = calcImageTransform(sceneGeometryIn,sceneGeometryOut,cameraOffsetPoint)
% Determine the rotation and translation matrices that describe the change
% in an image induced by the updated sceneParameters


pupilEllipseA1 = projectModelEye([ 0 1 0 1],sceneGeometryIn,'pupilRayFunc',[]);
pupilEllipseA2 = projectModelEye([-1 0 0 1],sceneGeometryIn,'pupilRayFunc',[]);
pupilEllipseA3 = projectModelEye([ 1 0 0 1],sceneGeometryIn,'pupilRayFunc',[]);
pupilEllipseB1 = projectModelEye([ 0 1 0 1],sceneGeometryOut,'pupilRayFunc',[]);
pupilEllipseB2 = projectModelEye([-1 0 0 1],sceneGeometryOut,'pupilRayFunc',[]);
pupilEllipseB3 = projectModelEye([ 1 0 0 1],sceneGeometryOut,'pupilRayFunc',[]);

A = [pupilEllipseA1(1:2)',pupilEllipseA2(1:2)',pupilEllipseA3(1:2)']-cameraOffsetPoint';
B = [pupilEllipseB1(1:2)',pupilEllipseB2(1:2)',pupilEllipseB3(1:2)']-cameraOffsetPoint';

regParams = absor(...
    A,...
    B,...
    'doScale',true,...
    'doTrans',true);
end


function displayImage = updateFrame(movingFrame,regParams,cameraOffsetPoint)

for ii=1:size(movingFrame,3)
    tmpImage = squeeze(movingFrame(:,:,ii));
    % Embed the movingFrame within a larger image that is padded
    % with mid-point background values
    padVals = round(size(tmpImage)./2);
    tmpImagePad = zeros(size(tmpImage)+padVals.*2)+125;
    tmpImagePad(padVals(1)+1:padVals(1)+size(tmpImage,1), ...
        padVals(2)+1:padVals(2)+size(tmpImage,2) ) = tmpImage;
    tmpImage = tmpImagePad;
    % Rotate the image
    tmpImage = imrotateAround(tmpImage, cameraOffsetPoint(2), cameraOffsetPoint(1), regParams.theta, 'bicubic');
    % Apply the x and y translation
    tmpImage = imtranslate(tmpImage,regParams.t','method','cubic');
    % Apply the scaling
    tmpImage = HardZoom(tmpImage,regParams.s);
    % Crop out the padding
    tmpImage = tmpImage(padVals(1)+1:padVals(1)+size(movingFrame,1), ...
        padVals(2)+1:padVals(2)+size(movingFrame,2));
    % Store this dimension
    displayImage(:,:,ii)=uint8(tmpImage);
end

end



function annotHandle = addAnnotation(text_str)
annotHandle = annotation('textbox',...
    [.80 .85 .1 .1],...
    'HorizontalAlignment','center',...
    'VerticalAlignment','middle',...
    'Margin',1,...
    'String',text_str,...
    'FontSize',9,...
    'FontName','Helvetica',...
    'EdgeColor',[1 1 1],...
    'LineWidth',1,...
    'BackgroundColor',[0.9  0.9 0.9],...
    'Color',[1 0 0]);
drawnow
end


function OutPicture = HardZoom(InPicture, ZoomFactor)
% Si el factor de escala es 1, no se hace nada
if ZoomFactor == 1
    OutPicture = InPicture;
    return;
end
% Se obtienen las dimensiones de las imágenes
ySize = size(InPicture, 1);
xSize = size(InPicture, 2);
zSize = size(InPicture, 3);
yCrop = floor(ySize / 2 * abs(ZoomFactor - 1));
xCrop = floor(xSize / 2 * abs(ZoomFactor - 1));
% Si el factor de escala es 0 se devuelve una imagen en negro
if ZoomFactor == 0
    OutPicture = uint8(zeros(ySize, xSize, zSize));
    return;
end
% Se reescala y se reposiciona en en centro
zoomPicture = imresize(InPicture, ZoomFactor);
ySizeZ = size(zoomPicture, 1);
xSizeZ = size(zoomPicture, 2);
if ZoomFactor > 1
    OutPicture = zoomPicture( 1+yCrop:yCrop+ySize, 1+xCrop:xCrop+xSize, :);
else
    OutPicture = uint8(zeros(ySize, xSize, zSize));
    OutPicture( 1+yCrop:yCrop+ySizeZ, 1+xCrop:xCrop+xSizeZ, :) = zoomPicture;
end
end


function output = imrotateAround(image, pointY, pointX, angle, varargin)
% ROTATEAROUND rotates an image.
%   ROTATED=ROTATEAROUND(IMAGE, POINTY, POINTX, ANGLE) rotates IMAGE around
%   the point [POINTY, POINTX] by ANGLE degrees. To rotate the image
%   clockwise, specify a negative value for ANGLE.
%
%   ROTATED=ROTATEAROUND(IMAGE, POINTY, POINTX, ANGLE, METHOD) rotates the
%   image with specified method:
%       'nearest'       Nearest-neighbor interpolation
%       'bilinear'      Bilinear interpolation
%       'bicubic'       Bicubic interpolation
%    The default is fast 'nearest'. Switch to 'bicubic' for nicer results.
%
%   Example
%   -------
%       imshow(rotateAround(imread('eight.tif'), 1, 1, 10));
%
%   See also IMROTATE, PADARRAY.
%   Contributed by Jan Motl (jan@motl.us)
%   $Revision: 1.2 $  $Date: 2014/05/01 12:08:01 $
% Parameter checking.
numvarargs = length(varargin);
if numvarargs > 1
    error('myfuns:somefun2Alt:TooManyInputs', ...
        'requires at most 1 optional input');
end
optargs = {'nearest'};    % Set defaults for optional inputs
optargs(1:numvarargs) = varargin;
[method] = optargs{:};    % Place optional args in memorable variable names
% Initialization.
[imageHeight, imageWidth, ~] = size(image);
centerX = floor(imageWidth/2+1);
centerY = floor(imageHeight/2+1);
dy = centerY-pointY;
dx = centerX-pointX;
% How much would the "rotate around" point shift if the
% image was rotated about the image center.
[theta, rho] = cart2pol(-dx,dy);
[newX, newY] = pol2cart(theta+angle*(pi/180), rho);
shiftX = round(pointX-(centerX+newX));
shiftY = round(pointY-(centerY-newY));
% Pad the image to preserve the whole image during the rotation.
padX = imageHeight;
padY = imageWidth;
padded = padarray(image, [padY padX],125);
% Rotate the image around the center.
rot = imrotate(padded, angle, method, 'crop');
% Crop the image.
output = rot(padY+1-shiftY:end-padY-shiftY, padX+1-shiftX:end-padX-shiftX, :);
end


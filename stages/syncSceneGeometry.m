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
%  'outputFileSuffix'     - Char vector. A string that is appended to the
%                           file name of the saved sceneGeometry and
%                           diagnostic plots.
%  'alignMethod'          - Cell or char vec. Controls how the routine
%                           selects frames from the pupilFile acquisition
%                           to use as the fixed target to which the
%                           sceneGeometry is adjusted. Valid options are:
%                               {'gazePre','gazePost','shape','gazeCalTest'}
%  'searchStrategy'       - Char vec. The searchStrategy is passed to the
%                           routine estimateSceneParams.m and defines which
%                           parameters and constraints ae used in the
%                           search. The default value of 'sceneSync' would
%                           only be changed in code and method development
%                           contexts.
%  'cameraTorsion','cameraDepth' - Scalar. If set, these values over-ride
%                           the x0 values that are taken from the
%                           sceneGeometry for the videoStemNameIn
%                           acquisition.
%  'frameSet'             - A 1xm vector that specifies the m indices
%                           (indexed from 1) that identify the set of
%                           frames from the acquisition to guide the search
%  'gazeTargets'          - A 2xm matrix that provides the positions, in
%                           degrees of visual angle, of fixation targets
%                           that correspond to each of the frames. Nan
%                           values are acceptable and indicate that the
%                           gaze position is unknown for a frame.
%
% Examples:
%{
%}


%% Parse vargin for options passed here
p = inputParser; p.KeepUnmatched = true; p.PartialMatching = false;

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

% Optional analysis params
p.addParameter('outputFileSuffix','',@ischar);
p.addParameter('alignMethod','gazePre',@(x)(ischar(x) | iscell(x)));
p.addParameter('searchStrategy','sceneSync',@ischar);
p.addParameter('usePriorResultAsX0',false,@islogical);
p.addParameter('cameraTorsion',[],@(x)(isempty(x) || isscalar(x)));
p.addParameter('cameraDepth',[],@(x)(isempty(x) || isscalar(x)));
p.addParameter('cameraTrans',[],@(x)(isempty(x) || isvector(x)));
p.addParameter('frameSet',[],@(x)(isempty(x) || isvector(x)));
p.addParameter('gazeTargets',[],@(x)(isempty(x) || ismatrix(x)));
p.addParameter('eyePoseLB',[-89,-89,0,0.1],@isnumeric);
p.addParameter('eyePoseUB',[89,89,0,5],@isnumeric);


%% Parse and check the parameters
p.parse(videoStemNameIn, videoStemNameOut, varargin{:});

% Handle early exit if either of the videoStems are empty
if isempty(videoStemNameIn) || isempty(videoStemNameOut)
    if p.Results.verbose
        fprintf('One or both of the videoStemNames was empty; returning.\n');
    end
    return
end

% Handle the alignMethod variable type
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

% Obtain information regarding the fixation frame
[fixIdxIn, ~, eyePoseFixationIn, rhoIn, thetaIn] = selectFrames.gazeCal(videoStemNameIn);

% Load in the video image for this frame.
videoFileNameIn = [videoStemNameIn '_gray.avi'];
if isfile(videoFileNameIn)
    videoFrameIn = makeMedianVideoImage(videoFileNameIn,'startFrame',fixIdxIn,'nFrames',1);
else
    videoFrameIn = [];
end

% Get the scene and eye parameters from the sceneGeometryIn
model.eye.x0 = sceneGeometryIn.meta.estimateSceneParams.xEye;
model.scene.x0 = sceneGeometryIn.meta.estimateSceneParams.xScene;

% Move any commonDepth value from the eye to the scene parameter set
depthIdx = find(strcmp(sceneGeometryIn.meta.estimateSceneParams.obj.model.scene.paramLabels,'depth'));
commonDepthIdx = find(strcmp(sceneGeometryIn.meta.estimateSceneParams.obj.model.eye.paramLabels,'commonDepth'));
model.scene.x0(depthIdx) = model.scene.x0(depthIdx) + model.eye.x0(commonDepthIdx);
model.eye.x0(commonDepthIdx) = 0;

% If cameraTorsion, cameraDepth, or cameraTrans has been passed, use this
% to over-write the default value in the model parameters
if ~isempty(p.Results.cameraTorsion)
    model.scene.x0(strcmp(sceneGeometryIn.meta.estimateSceneParams.obj.model.scene.paramLabels,'torsion')) = p.Results.cameraTorsion;
end
if ~isempty(p.Results.cameraDepth)
    model.scene.x0(strcmp(sceneGeometryIn.meta.estimateSceneParams.obj.model.scene.paramLabels,'depth')) = p.Results.cameraDepth;
end
if ~isempty(p.Results.cameraTrans)
    model.scene.x0(strcmp(sceneGeometryIn.meta.estimateSceneParams.obj.model.scene.paramLabels,'horiz')) = p.Results.cameraTrans(1);
    model.scene.x0(strcmp(sceneGeometryIn.meta.estimateSceneParams.obj.model.scene.paramLabels,'vert')) = p.Results.cameraTrans(2);
end

% Define the filename for the sceneGeometry output
sceneGeometryNameOut = [videoStemNameOut '_sceneGeometry' p.Results.outputFileSuffix '.mat'];

% If usePriorResultAsX0 is set to true, and there is a prior result, load
% the prior sceneGeometry result and set the model x0 values with the prior
% results
if p.Results.usePriorResultAsX0 && exist(sceneGeometryNameOut, 'file') == 2
    % Load the sceneGeometryOut variable into memory
    dataLoad=load(sceneGeometryNameOut);
    sceneGeometryOut=dataLoad.sceneGeometry;
    clear dataLoad
    % Check that there is a prior estimateSceneParams field
    if isfield(sceneGeometryOut.meta,'estimateSceneParams')
        model.scene.x0 = sceneGeometryOut.meta.estimateSceneParams.xScene;
        model.head.x0 = sceneGeometryOut.meta.estimateSceneParams.xHead;
    end
    clear sceneGeometryOut
end


% Obtain the eye and error args
eyeArgs = sceneGeometryIn.meta.estimateSceneParams.obj.setupArgs;
errorArgs = { ...
    'poseRegParams',sceneGeometryIn.meta.estimateSceneParams.obj.poseRegParams, ...
    'eyePoseLB',p.Results.eyePoseLB,'eyePoseUB',p.Results.eyePoseUB};


%% Select frames to guide the search
% We need a set of frames that have a high-quality pupil perimeter and an
% associated glint location. These frames should be distributed across
% time, and capture the eye in a variety of poses. The sceneGeometry will
% be adjusted to best fit these frames. Also, we would ideally have a frame
% in which the eye is posed such that it is looking at the center point of
% a fixation array (gazeTarget [0;0]). If we can find such a frame, it will
% be the first entry in the list.

if ~isempty(p.Results.frameSet)
    
    % A frameSet was passed; use that
    frameSet = p.Results.frameSet;
    
    % Use the gazeTargets if passed, otherwise create a matrix of nans
    if ~isempty(p.Results.gazeTargets)
        gazeTargets = p.Results.gazeTargets;
    else
        gazeTargets = nan(2,length(frameSet));
    end
    
    % Because the frameSet was defined, the alignMethod is not used
    alignMethod = 'passedSet';

else
    % We will create a frameSet
    nFramesToReturn = 15;
    
    switch alignMethod
        case 'gazePre'
            
            % This frame is for display and fixation
            [frameSet, gazeTargets] = selectFrames.gazePre(videoStemNameOut);
            
        case 'gazePost'
            
            % This frame is for display and fixation
            [frameSet, gazeTargets] = selectFrames.gazePost(videoStemNameOut);
            
        case 'shape'
            
            % This frame is for display only
            [frameSet, gazeTargets] = selectFrames.shape(videoStemNameOut, rhoIn, thetaIn);
            
        case 'gazeCalTest'
            
            % Get all of the gaze target frames from the source. They will
            % be ordered such that the fixation ([0;0]) gaze position is
            % first.
            [frameSet, gazeTargets] = selectFrames.gazeCalTest(videoStemNameOut);
            
            % We need fewer additional frames, as we already have 9
            nFramesToReturn = 10;
            
        otherwise
            error('This is not a defined align method')
    end
    
    % Add frames that are distributed across time and space. Start with a
    % low distValThresh and increase as needed to reach the desired frames
    % or until the maximum suitable threshold is reached.
    stillSearching = true;
    distValsThreshold = 0.2;
    maxDistValsThreshold = 0.35;
    while stillSearching
        [frameSetA, gazeTargetsA] = selectFrames.gridTime(videoStemNameOut,'nFramesToReturn',nFramesToReturn,'distValsThreshold',distValsThreshold);
        [frameSetB, gazeTargetsB] = selectFrames.gridSpace(videoStemNameOut,'nFramesToReturn',nFramesToReturn,'distValsThreshold',distValsThreshold);
        if length(frameSetA)>=nFramesToReturn && length(frameSetB)>=nFramesToReturn
            stillSearching = false;
        else
            distValsThreshold = distValsThreshold + 0.025;
        end
        if distValsThreshold > maxDistValsThreshold
            stillSearching = false;
        end
    end
    
    % Assemble the frames
    frameSet = [frameSet frameSetA frameSetB];
    gazeTargets = [gazeTargets gazeTargetsA gazeTargetsB];
    
    % Remove duplicate frames
    [~, uniqueFrames] = unique(frameSet,'stable');
    frameSet = frameSet(uniqueFrames);
    gazeTargets = gazeTargets(:,uniqueFrames);
    
end


% Load in the video image for this frame.
videoFrameOut = makeMedianVideoImage([videoStemNameOut '_gray.avi'],'startFrame',frameSet(1),'nFrames',1);

% If the VideoFrameIn is empty, create a blank frame
if isempty(videoFrameIn)
    videoFrameIn = zeros(size(videoFrameOut));
end

%% Perform the synchronization search
estimateSceneParams(videoStemNameOut, frameSet, gazeTargets, ...
    'outputFileSuffix',p.Results.outputFileSuffix,...
    'searchStrategy',p.Results.searchStrategy,...
    'model',model,...
    'eyeArgs',eyeArgs,'errorArgs',errorArgs, ...
    'verbose',p.Results.verbose,...
    'useParallel',p.Results.useParallel,...
    'nWorkers',p.Results.nWorkers);

% Load the sceneGeometryOut variable into memory
dataLoad=load(sceneGeometryNameOut);
sceneGeometryOut=dataLoad.sceneGeometry;
clear dataLoad

% Store the eyePose for the fixation frame of sceneGeometryOut
eyePoseFixationOut = sceneGeometryOut.meta.estimateSceneParams.obj.modelEyePose(1,:);

% Obtain the camera translation for the fixation frame of sceneGeometryOut
cameraTrans = sceneGeometryOut.meta.estimateSceneParams.obj.modelCameraTrans(:,1);

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
    str = strsplit(videoStemNameIn,filesep);
    text(20,30,str{end}, 'Color', 'g','Fontsize',16,'Interpreter','none');
    msg = ['frame ' num2str(fixIdxIn)];
    addAnnotation(msg);
    hold on
    plot([size(displayImage,2)/2, size(displayImage,2)/2],[0 size(displayImage,1)],'-b');
    plot([0 size(displayImage,1)],[size(displayImage,1)/2, size(displayImage,1)/2],'-b');
    tmpFrame = getframe(gcf);
    imageSet(1) = {tmpFrame.cdata};
    close(tmpFig);
    
    
    %% videoStemNameOut
    % The acquisition for which we have new sceneGeometry
    
    displayImage = videoFrameOut;
    tmpFig = figure('visible','off');
    renderEyePose(eyePoseFixationOut, sceneGeometryOut, ...
        'cameraTrans',cameraTrans, ...
        'newFigure', false, 'visible', false, ...
        'backgroundImage',displayImage, ...
        'showAzimuthPlane', true, ...
        'modelEyeLabelNames', {'retina' 'pupilEllipse' 'cornea' 'glint_01'}, ...
        'modelEyePlotColors', {'.w' '-g' '.y' 'xr'}, ...
        'modelEyeSymbolSizeScaler',1.5,...
        'modelEyeAlpha', [0.25 0.25 0.25 1]);
    str = strsplit(videoStemNameOut,filesep);
    text(20,30,str{end}, 'Color', 'r','Fontsize',16,'Interpreter','none');
    msg = ['frame ' num2str(frameSet(1))];
    addAnnotation(msg);
    % Add cross hairs
    hold on
    plot([size(displayImage,2)/2, size(displayImage,2)/2],[0 size(displayImage,1)],'-b');
    plot([0 size(displayImage,1)],[size(displayImage,1)/2, size(displayImage,1)/2],'-b');
    tmpFrame = getframe(gcf);
    imageSet(2) = {tmpFrame.cdata};
    close(tmpFig);
    
    
    %% Difference image
    
    % Transform the videoFrameIn to the space of videoFrameOut
    regParams = calcImageTransform(sceneGeometryIn,sceneGeometryOut,cameraTrans,cameraOffsetPoint);
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
    
    % Add a text summary below. If any delta fixation angle is greater than
    % 1 deg, and the align method was not "shape", print the message text
    % in red to alert that this was a large eye rotation change.
    deltaX = sceneGeometryOut.meta.estimateSceneParams.xScene - model.scene.x0;
    deltaPose = eyePoseFixationOut - eyePoseFixationIn;
    
    msg = sprintf('delta torsion [deg] = %2.1f',deltaX(3));
    annotation('textbox', [0.5, .175, 0, 0], 'string', msg,'FitBoxToText','on','LineStyle','none','HorizontalAlignment','center','Interpreter','none')
    msg = sprintf('delta translation [mm] [x; y; z] = [%2.3f; %2.3f; %2.3f]',deltaX(4:6));
    annotation('textbox', [0.5, .125, 0, 0], 'string', msg,'FitBoxToText','on','LineStyle','none','HorizontalAlignment','center','Interpreter','none')
    msg = sprintf('delta eye pose [azi, ele, tor, radius] = [%2.3f, %2.3f, %2.3f, %2.3f]',deltaPose);
    msgColor = 'black';
    if any(abs(deltaPose(1:3)) > 1.0) && ~strcmp(alignMethod,'shape')
        msgColor = 'red';
    end
    annotation('textbox', [0.5, .075, 0, 0], 'string', msg,'Color',msgColor,'FitBoxToText','on','LineStyle','none','HorizontalAlignment','center','Interpreter','none')
    
    % Save and close the figure
    print(figHandle,[videoStemNameOut '_sceneGeometry_sceneSync_QA' p.Results.outputFileSuffix '.pdf'],'-dpdf','-bestfit');
    close(figHandle);
    
end



end % Main function





%% LOCAL FUNCTIONS


function regParams = calcImageTransform(sceneGeometryIn,sceneGeometryOut,cameraTrans,cameraOffsetPoint)
% Determine the rotation and translation matrices that describe the change
% in an image induced by the updated sceneParameters

% Adjust the sceneGeometryOut to account for the camera translation
sceneGeometryOut.cameraPosition.translation = ...
    sceneGeometryOut.cameraPosition.translation + cameraTrans;

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
% Se obtienen las dimensiones de las im�genes
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


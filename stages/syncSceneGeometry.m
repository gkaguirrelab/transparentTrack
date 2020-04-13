function syncSceneGeometry(videoStemNameIn, videoStemNameOut, varargin)
% Update sceneGeometry camera position to match an acquisition
%
% Syntax:
%  syncSceneGeometry(videoStemNameIn, videoStemNameOut, varargin)
%
% Description:
%   A sceneGeometry file is created for a given acquisition. Included in
%   the sceneGeometry is a specification of properties of the extrinsic
%   camera matrix, including the position of the camera in space relative
%   to the coordinates, which have as their origin the anterior surface of
%   the cornea along the optical axis of the eye. If we wish to use this
%   sceneGeometry file for the analysis of data from other acqusitions for
%   a given subject, we need to deal with the possibility that the subject
%   has moved their head between acquisitions. As the coordinate system is
%   based upon a fixed anatomical landmark of the eye, the effect of head
%   translation in this system is to change the camera position. This
%   routine assists in calculating an updated camera position for a given
%   acquisition.
%
%   Select a sceneGeometry with the UI file picker, then select the
%   acquisitions to which you would like to align the sceneGeometry.
%   An image derived from the sceneGeometry and the acquisition is shown,
%   and controls are used to adjust the acquisition to match the
%   sceneGeometry.
%
% Inputs:
%   pupilFileName         - Full path to the .mat file that contains the
%                           pupil data to which the sceneGeometry should be
%                           synced.
%   
% Optional key/value pairs (display and I/O):
%  'sceneGeometryFileNameToSync' - Full path to the .mat file that contains
%                           the sceneGeometry to be used.
%  'displayMode'          - Logical. Controls if a GUI interface is
%                           provided after the search stage to allow the
%                           user to adjust the scene geometry parameters
%                           by hand.
%  'verbose'              - Logical.
%  'saveAdjustedSceneGeometry' - Logical. Controls if the adjusted 
%                           sceneGeometry file is saved. If so, the file
%                           name will be derived from the pupilFile name.
%  'saveDiagnosticPlot'   - Logical.
%  'doNotSyncSceneToItself' - Logical. Strange things could happen if a the
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
%  'alignMethod'          - Char vec. Controls how the routine selects
%                           frames from the pupilFile acquisition to use
%                           as the fixed target to which the sceneGeometry
%                           is adjusted. Valid options are:
%                               {'gazePre','gazePost','shape'}
%  'sceneSyncX'           - 1x8 vector. If set, these scene parameter
%                           values will be used for the adjusted
%                           sceneGeometry, instead of performing a search
%                           to try and find optimal parameters.
%  'sceneSyncBound'       - 1x8 vector. Defines the +-bounds on the scene
%                           parameters for the search that brings the
%                           source sceneGeometry into alignment with the
%                           fixed, target pupilFile data. As properties of
%                           the eye itself (rotation center, corneal
%                           curvature) are considered fixed, the last four
%                           values of this vector are usually set to zero.
%                           Also, given a single eyePose as a target, the
%                           model has difficulty adjusting torsion (the
%                           first parameter) and depth (the fourth
%                           parameter) so these are typically set to zero
%                           as well.
%  'eyePositionTargetLengthFrames' - Scalar. The number of sequential 
%                           frames from the target acquisition that will be
%                           found and used to define the position of the
%                           eye to be fit.
%  'gazeErrorThreshTol'   - Scalar. The run of frames must have a deviation
%                           of less than this value. The precise meaning of
%                           the value will differ for the different
%                           alignment methods.
%
% Examples:
%{
    % Invoke the file picker GUI
    syncSceneGeometry('','displayMode',true,'verbose',true,'alignMethod','shape');
%}
%{
    sceneGeometryFileNameToSync = '/Users/aguirre/Dropbox (Aguirre-Brainard Lab)/TOME_processing/session2_spatialStimuli/TOME_3015/032417/EyeTracking/GazeCal03_sceneGeometry.mat';
    pupilFileName = '/Users/aguirre/Dropbox (Aguirre-Brainard Lab)/TOME_processing/session2_spatialStimuli/TOME_3015/032417/EyeTracking/tfMRI_MOVIE_AP_run01_pupil.mat';
    syncSceneGeometry(pupilFileName, 'sceneGeometryFileNameToSync', sceneGeometryFileNameToSync) 
%}

%% Parse vargin for options passed here
p = inputParser; p.KeepUnmatched = true;

% Required
p.addRequired('videoStemNameIn',@ischar);
p.addRequired('videoStemNameOut',@ischar);

% Optional display and I/O params
p.addParameter('displayMode',false,@islogical);
p.addParameter('verbose',true,@islogical);
p.addParameter('saveAdjustedSceneGeometry',true,@islogical);
p.addParameter('saveDiagnosticPlot',true,@islogical);
p.addParameter('doNotSyncSceneToItself',true,@islogical);

% Optional flow control params
p.addParameter('useParallel',true,@islogical);
p.addParameter('nWorkers',[],@(x)(isempty(x) || isnumeric(x)));

% Optional fitting params
p.addParameter('alignMethod','gazePre',@(x)(ischar(x) | iscell(x)));


%% Parse and check the parameters
p.parse(videoStemNameIn, videoStemNameOut, varargin{:});


%% Check if we are syncing a sceneGeometry to itself
if p.Results.doNotSyncSceneToItself && strcmp(videoStemNameIn,videoStemNameOut)
    if p.Results.verbose
        fprintf('Detected that the sceneGeometry source is the same as the target acquisition; returning.\n');
    end
    return
end


%% Announce we are starting
ticObject = tic();
if p.Results.verbose
    fprintf(['Syncing scene geometry. Started ' char(datetime('now')) '\n']);
end


%% Load sceneGeometry for the videoStemNameIn
% This sceneGeometry provides the eye model and the parameters that
% transform eyePose into gaze position.

% Load the sceneGeometry variable into memory
dataLoad=load([videoStemNameIn '_sceneGeometry.mat');
sceneGeometryIn=dataLoad.sceneGeometry;
clear dataLoad

% Get the camera offset point. This is where the center of the video image
% lands on the camera sensor array. We will need this later when we make
% some diagnostic images of the alignment of the source and target videos.
cameraOffsetPoint = [sceneGeometryIn.cameraIntrinsic.matrix(1,3), ...
    sceneGeometryIn.cameraIntrinsic.matrix(2,3)];


%% Find the fixation frame for the source image
% This is the sceneGeometry file that (typically) was derived during a
% gazeCalibration procedure.

% Which of the list of frames is the [0,0] fixation frame
idx = find((sceneGeometryIn.meta.estimateSceneParams.gazeTargets(1,:)==0).*(sceneGeometryIn.meta.estimateSceneParams.gazeTargets(2,:)==0));

% Store the eyePose for this frame
eyePoseFixationIn = sceneGeometryIn.meta.estimateSceneParams.obj.modelEyePose(idx,:);

% Store the pupilEllipse for this frame
pupilEllipseFixationIn = sceneGeometryIn.meta.estimateSceneParams.obj.modelPupilEllipse(idx,:);

% Find the shape of the pupil for this frame, expressed as theta and rho
% values (SEE: csaEllipseError)
rhoIn = 1-sqrt(1-pupilEllipseFixationIn(4)^2);
thetaIn = pupilEllipseFixationIn(5)*2;

% Get the scene model parameters
model.eye.x0 = sceneGeometryIn.meta.estimateSceneParams.xEye;
eyeArgs = sceneGeometryIn.meta.estimateSceneParams.obj.setupArgs;
errorArgs = { ...
        'poseRegParams',sceneGeometryIn.meta.estimateSceneParams.obj.poseRegParams,...
        'vecRegParams',sceneGeometryIn.meta.estimateSceneParams.obj.vectorRegParams};

% Get the cameraDepth from the source sceneGeometry file
cameraDepth = sceneGeometryIn.meta.estimateSceneParams.xScene(end);
    
% Load in the video image for this frame.
tmp = fullfile(sceneGeometryInPath,[videoStemNameIn '_gray.avi']);
absIdx = sceneGeometryIn.meta.estimateSceneParams.frameSet(idx);
videoFrameIn = makeMedianVideoImage(tmp,'startFrame',absIdx,'nFrames',1);


%% Identify frames from videoStemNameOut to guide the search
% Every search includes frames that are distributed across time and space.
% Additionally, we wish to identify a frame in the videoStemNameOut
% acquisition that best corresponds to the eye when it is looking at the
% fixation point.

[frameSetA, gazeTargetsA] = selectFrames.gridTime(videoStemNameOut);
[frameSetB, gazeTargetsB] = selectFrames.gridSpace(videoStemNameOut);

switch p.Results.alignMethod
    case 'gazePre'
        [frameSetC, gazeTargetsC] = selectFrames.gazePre(videoStemNameOut);
    case 'gazePost'
        [frameSetC, gazeTargetsC] = selectFrames.gazePost(videoStemNameOut);
    case 'shape'
        [frameSetC, gazeTargetsC] = selectFrames.shape(videoStemNameOut, rhoIn, thetaIn);    
end

frameSet = [frameSetA frameSetB frameSetC];
gazeTargets = [gazeTargetsA gazeTargetsB gazeTargetsC];


%% Perform the synchronization search
estimateSceneParams(videoStemNameOut, frameSet, gazeTargets, ...
    'searchStrategy','gazeCalReg','cameraDepth',cameraDepth,'model',model,...
    'eyeArgs',eyeArgs,'errorArgs',errorArgs);

% 
% 
% %% Optional display mode
% % If we are in display mode, offer an interface for the user to manually
% % adjust the delta variables
% 
% if p.Results.displayMode
%         
%     % Provide some instructions for the operator
%     fprintf([sceneGeometryInPath '\n']);
%     fprintf('Adjust horizontal and vertical camera translation with the arrow keys.\n');
%     fprintf('Adjust depth camera translation with + and -.\n');
%     fprintf('Adjust camera torsion with j and k.\n');
%     fprintf('Switch between moving and fixed image by pressing f.\n');
%     fprintf('Turn on and off perimeter display with p.\n');
%     fprintf('Turn on and off model display with m.\n');
%     fprintf('Press esc to exit.\n\n');
%     
%     % Create a figure
%     figHandle = figure();
%     imshow(videoFrameFixed,[],'Border','tight');
%     ax = gca;
%     ax.Toolbar = [];
%     hold on
%     text(20,30,'FIXED', 'Color', 'g','Fontsize',16);
%     
%     % Prepare for the loop
%     showMoving = false;
%     showPerimeter=false;
%     showModel=false;
%     stillWorking = true;
%     
%     % Enter the while stillWorking loop
%     while stillWorking
%         
%         % Prepare to update the image
%         hold off
%         
%         if showMoving
%             % Work with the moving frame
%             displayImage = videoFrameMoving;
%             
%             % Update the perimeter points
%             XpDisplay = XpMoving; YpDisplay = YpMoving;
%                         
%         else
%             % Work with the fixed frame
%             displayImage = videoFrameFixed;
%             
%             % Update the perimeter
%             XpDisplay = XpFixed; YpDisplay = YpFixed;
%             
%         end
%         
%         % Display the perimeter points
%         if showPerimeter
%             idx = sub2ind(size(displayImage),round(YpDisplay),round(XpDisplay));
%             displayImage(idx)=255;
%         end
%         
%         % Add the eye model
%         if showModel
%             % Let the user know this will take a few seconds
%             text_str = 'Updating model...';
%             annotHandle = addAnnotation(text_str);
%             % Set the eye pose, depending upon if we are looking at the
%             % fixed or moving image
%             if showMoving
%                 sceneGeometryDisplay = sceneGeometryIn;
%                 eyePoseDisplay = eyePoseFixationIn;
%             else
%                 sceneGeometryDisplay = updateSceneGeometry( sceneGeometryIn, x );
%                 eyePoseDisplay = eyePoseEllipseFit(XpFixed, YpFixed, sceneGeometryDisplay,'glintCoord',glintCoordFixed,keyVals{:});
%             end
%             % Render the eye model
%             [tmpHandle,~,displayImage]=renderEyePose(eyePoseDisplay, sceneGeometryDisplay, ...
%                 'newFigure', true, 'visible', false, ...
%                 'backgroundImage', displayImage, ...
%                 'showAzimuthPlane', true, ...
%                 'modelEyeLabelNames', {'retina' 'pupilEllipse' 'cornea' 'glint_01'}, ...
%                 'modelEyePlotColors', {'.w' '-g' '.y' 'xr'}, ...
%                 'modelEyeSymbolSizeScaler',1.5,...
%                 'modelEyeAlpha', [0.25 0.25 0.25 1]);
%             close(tmpHandle);
%             displayImage = displayImage.cdata;
%             % Remove the updating annotation
%             delete(annotHandle);
%         end
% 
%         % Move the moving image
%         if showMoving
%             % Get the 2D registration params that brings the movingImage
%             % into the fixedImage space
%             regParams = calcImageTransform(sceneGeometryIn,x,cameraOffsetPoint);
%             
%             % Update the image            
%             displayImage = updateFrame(displayImage,regParams,cameraOffsetPoint);
%         end
%             
%         % Display the image
%         imshow(displayImage,[],'Border','tight');
%         ax = gca;
%         ax.Toolbar = [];
%         hold on
%         % Add a label
%         if showMoving
%             text(20,30,'MOVING', 'Color', 'r','Fontsize',16);
%         else
%             text(20,30,'FIXED', 'Color', 'g','Fontsize',16);
%         end
%         
%         % Add a marker for the camera CoP
%         plot(cameraOffsetPoint(1),cameraOffsetPoint(2),'+c');
%         
%         keyAction = waitforbuttonpress;
%         if keyAction
%             keyChoiceValue = double(get(gcf,'CurrentCharacter'));
%             switch keyChoiceValue
%                 case 28
%                     x(2)=x(2)-0.1;
%                 case 29
%                     x(2)=x(2)+0.1;
%                 case 30
%                     x(3)=x(3)-0.1;
%                 case 31
%                     x(3)=x(3)+0.1;
%                 case {45 95}
%                     x(4)=x(4)-1;
%                 case {61 43}
%                     x(4)=x(4)+1;
%                 case 102
%                     showMoving = ~showMoving;
%                 case 106
%                     x(1)=x(1)+0.5;
%                 case 107
%                     x(1)=x(1)-0.5;
%                 case 112
%                     showPerimeter = ~showPerimeter;
%                 case 109
%                     showModel = ~showModel;
%                 case 27
%                     text_str = 'finishing...';
%                     addAnnotation(text_str);
%                     stillWorking = false;
%                 otherwise
%                     % Nothing to do
%             end
%             
%         end
%     end
%     
%     close(figHandle);
% end
% 
% 
% %% Create the adjusted sceneGeometry
% sceneGeometryAdjusted = updateSceneGeometry( sceneGeometryIn, x );
% 
% 
% %% Update the fixation angles
% [~,modeledEyePose] = calcGlintGazeError( sceneGeometryAdjusted, args{:}, keyVals{:} );
% medianEyePoseFixed = median(modeledEyePose);
% sceneGeometryAdjusted.screenPosition.fixationEyePose = medianEyePoseFixed(1:2)';
% 
% 
% %% Create and save a diagnostic figure
% if p.Results.saveDiagnosticPlot
%     
%     % Moving frame -- Typically the gazeCal source
%     displayImage = videoFrameMoving;
%     tmpFig = figure('visible','off');
%     renderEyePose(eyePoseFixationIn, sceneGeometryIn, ...
%         'newFigure', false, 'visible', false, ...
%         'backgroundImage',displayImage, ...
%         'showAzimuthPlane', true, ...
%                 'modelEyeLabelNames', {'retina' 'pupilEllipse' 'cornea' 'glint_01'}, ...
%                 'modelEyePlotColors', {'.w' '-g' '.y' 'xr'}, ...
%                 'modelEyeSymbolSizeScaler',1.5,...
%                 'modelEyeAlpha', [0.25 0.25 0.25 1]);            
%     text(20,30,videoStemNameIn, 'Color', 'g','Fontsize',16,'Interpreter','none');
%     msg = ['frame ' num2str(absIdx)];
%     addAnnotation(msg);
%     hold on
%     plot([size(displayImage,2)/2, size(displayImage,2)/2],[0 size(displayImage,2)],'-b');
%     plot([0 size(displayImage,1)],[size(displayImage,1)/2, size(displayImage,1)/2],'-b');
%     tmpFrame = getframe(gcf);
%     imageSet(1) = {tmpFrame.cdata};
%     close(tmpFig);
% 
%         
%     % Fixed frame -- The acquisition for which we have new sceneGeometry
%     displayImage = videoFrameFixed;
%     eyePoseFixed = eyePoseEllipseFit(XpFixed, YpFixed, sceneGeometryAdjusted,'glintCoord',glintCoordFixed,keyVals{:});
%     tmpFig = figure('visible','off');
%     renderEyePose(eyePoseFixed, sceneGeometryAdjusted, ...
%         'newFigure', false, 'visible', false, ...
%         'backgroundImage',displayImage, ...
%         'showAzimuthPlane', true, ...
%                 'modelEyeLabelNames', {'retina' 'pupilEllipse' 'cornea' 'glint_01'}, ...
%                 'modelEyePlotColors', {'.w' '-g' '.y' 'xr'}, ...
%                 'modelEyeSymbolSizeScaler',1.5,...
%                 'modelEyeAlpha', [0.25 0.25 0.25 1]);            
%     text(20,30,sceneGeometryOutStem, 'Color', 'r','Fontsize',16,'Interpreter','none');
%     msg = ['frame ' num2str(referenceFrameFixed)];
%     addAnnotation(msg);
%     % Add cross hairs
%     hold on
%     plot([size(displayImage,2)/2, size(displayImage,2)/2],[0 size(displayImage,2)],'-b');
%     plot([0 size(displayImage,1)],[size(displayImage,1)/2, size(displayImage,1)/2],'-b');
%     tmpFrame = getframe(gcf);
%     imageSet(2) = {tmpFrame.cdata};
%     close(tmpFig);
%     
%     
%     % Difference image
%     regParams = calcImageTransform(sceneGeometryIn,x,cameraOffsetPoint);
%     adjMovingFrame = updateFrame(videoFrameMoving,regParams,cameraOffsetPoint);
%     displayImage = videoFrameFixed - double(adjMovingFrame);
%     tmpFig = figure('visible','off');
%     imshow(displayImage,[], 'Border', 'tight');
%     text(20,30,'Difference', 'Color', 'w','Fontsize',16,'Interpreter','none');
%     tmpFrame = getframe(gcf);
%     imageSet(3) = {tmpFrame.cdata};
%     close(tmpFig);
%     
%     % Prepare the figure
%     figHandle=figure('visible','off');
%     set(gcf,'PaperOrientation','portrait');
%     
%     set(figHandle, 'Units','inches')
%     height = 4;
%     width = 12;
%     
%     % The last two parameters of 'Position' define the figure size
%     set(figHandle, 'Position',[25 5 width height],...
%         'PaperSize',[width height],...
%         'PaperPositionMode','auto',...
%         'Color','w',...
%         'Renderer','painters'...
%         );
%     
%     % Post the montage of the imageSet
%     montage(imageSet,'Size', [1 3]);
%     
%     % Post the title
%     pathParts = strsplit(sceneGeometryInPath,filesep);
%     titleString = [fullfile(pathParts{end-4:end-2})];
%     title(titleString,'Interpreter','none')
%     
%     % Report the alignment method
%     annotation('textbox', [0.15, .125, 0, 0], 'string', alignMethod,'FontWeight','bold','FitBoxToText','on','LineStyle','none','HorizontalAlignment','left','Interpreter','none')     
%     
%     % Report the run lengths
%     msg = sprintf('runLength fixed = %2.0f',runLengthFixed);
%     annotation('textbox', [0.75, .125, 0, 0], 'string', msg,'FitBoxToText','on','LineStyle','none','HorizontalAlignment','left','Interpreter','none')     
%     
%     % Add a text summary below. If any delta fixation angle is geater than
%     % 1 deg, print the message text in red to alert that this was a large
%     % eye rotation change.
%     deltaX = x-x0;
%     deltaPose = medianEyePoseFixed - eyePoseFixationIn;
%     msg = sprintf('delta torsion [deg] = %2.1f',deltaX(1));
%     annotation('textbox', [0.5, .175, 0, 0], 'string', msg,'FitBoxToText','on','LineStyle','none','HorizontalAlignment','center','Interpreter','none')
%     msg = sprintf('delta translation [mm] [x; y; z] = [%2.3f; %2.3f; %2.3f]',deltaX(2:4));
%     annotation('textbox', [0.5, .125, 0, 0], 'string', msg,'FitBoxToText','on','LineStyle','none','HorizontalAlignment','center','Interpreter','none')
%     msg = sprintf('delta eye pose [azi, ele, tor, radius] = [%2.3f, %2.3f, %2.3f, %2.3f]',deltaPose);
%     msgColor = 'black';
%     if any(abs(deltaPose(1:3)) > 0) 
%         msgColor = 'red';
%     end
%     annotation('textbox', [0.5, .075, 0, 0], 'string', msg,'Color',msgColor,'FitBoxToText','on','LineStyle','none','HorizontalAlignment','center','Interpreter','none')
%     
%     % Save and close the figure
%     tmp = fullfile(sceneGeometryOutPath,[sceneGeometryOutStem '_sceneSync_QA.pdf']);
%     print(figHandle,tmp,'-dpdf','-bestfit');
%     close(figHandle);
%     
% end
% 
% % Get the execution time
% executionTime = toc(ticObject);
% 
% %% Save the adjusted sceneGeometry
% if p.Results.saveAdjustedSceneGeometry
%     
%     % If the movingFrame is from a point after the start of the
%     % acquisition, adjust the camera translation position once more to
%     % account for any head motion that took place between the start of the
%     % acquisition and the reference frame
%     sceneGeometryAdjusted.cameraPosition.translation = ...
%         sceneGeometryAdjusted.cameraPosition.translation + ...
%         relativeCameraPosition.values(:,referenceFrameFixed);
%     
%     % Add the meta data
%     sceneGeometryAdjusted.meta.syncSceneGeometry = p.Results;
%     sceneGeometryAdjusted.meta.syncSceneGeometry.x = x;
%     sceneGeometryAdjusted.meta.executionTime = executionTime;
%     
%     % Set the variable name
%     sceneGeometry = sceneGeometryAdjusted;
%     
%     % Save it
%     tmp = fullfile(sceneGeometryOutPath,[sceneGeometryOutStem '_sceneGeometry.mat']);
%     save(tmp,'sceneGeometry');
% end
% 
% % If we are in display mode, report the values to the console
% if p.Results.displayMode
%     tmp=strsplit(sceneGeometryInPath,filesep);
%     outline = [tmp{end-3} char(9) tmp{end-2} char(9) 'fixed: ' videoStemNameIn ', moving: ' sceneGeometryOutStem '\n'] ;
%     fprintf(outline)
%     outline = ['alignMethod' char(9) 'x' char(9) 'eyePose\n'];
%     fprintf(outline)
%     outline = sprintf(['{''' alignMethod '''}' char(9) '[ %2.2f, %2.2f, %2.2f, %2.2f, %2.2f, %2.2f, %2.2f, %2.2f ]' char(9) '[ %2.2f, %2.2f, %2.2f, %2.2f ]\n'],x,medianEyePoseFixed);
%     fprintf(outline)
%     fprintf('\n')
% end
% 
% 
% %% alert the user that we are done with the routine
% if p.Results.verbose
%     executionTime
%     fprintf('\n');
% end


end % Main function





%% LOCAL FUNCTIONS


function regParams = calcImageTransform(sceneGeometry,x,cameraOffsetPoint)
% Determine the rotation and translation matrices that describe the change
% in an image induced by the updated sceneParameters
f = updateSceneGeometry(sceneGeometry,x);

pupilEllipseA1 = projectModelEye([ 0 1 0 1],sceneGeometry,'pupilRayFunc',[]);
pupilEllipseA2 = projectModelEye([-1 0 0 1],sceneGeometry,'pupilRayFunc',[]);
pupilEllipseA3 = projectModelEye([ 1 0 0 1],sceneGeometry,'pupilRayFunc',[]);
pupilEllipseB1 = projectModelEye([ 0 1 0 1],f,'pupilRayFunc',[]);
pupilEllipseB2 = projectModelEye([-1 0 0 1],f,'pupilRayFunc',[]);
pupilEllipseB3 = projectModelEye([ 1 0 0 1],f,'pupilRayFunc',[]);

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


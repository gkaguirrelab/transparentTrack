function syncSceneGeometry(sceneGeometryFileName, pupilFileName, varargin)
% The change in camera position between a sceneGeometry and an acquisition
%
% Syntax:
%  
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
%   sceneGeometryFileName - Full path to the .mat file that contains the
%                           sceneGeometry to be used.
%   pupilFileName         - Full path to the .mat file that contains the
%                           pupil data to which the sceneGeometry should be
%                           synced.
%
% Examples:
%{
    % Invoke the file picker GUI
    syncSceneGeometry('','');
%}


%% Parse vargin for options passed here
p = inputParser; p.KeepUnmatched = true;

% Required
p.addRequired('sceneGeometryFileName',@ischar);
p.addRequired('pupilFileName',@ischar);

% Optional display and I/O params
p.addParameter('verbose',false,@islogical);
p.addParameter('displayMode',false,@islogical);
p.addParameter('saveAdjustedSceneGeome',true,@islogical);

% Optional fitting params
p.addParameter('alignMethod','gaze',@ischar);
p.addParameter('adjustedCameraPositionTranslation',[],@isnumeric);
p.addParameter('adjustedCameraPositionTorsion',[],@isnumeric);
p.addParameter('adjustedCameraPositionFixationAngles',[],@isnumeric);




%% Parse and check the parameters
p.parse(sceneGeometryFileName, pupilFileName, varargin{:});



%% Load the sceneGeometry file
if isempty(sceneGeometryFileName)
    % Open a file picker UI to select a sceneGeometry
    [file,path] = uigetfile(fullfile('.','*_sceneGeometry.mat'),'Choose a sceneGeometry file');
    sceneGeometryIn = fullfile(path,file);
else
    sceneGeometryIn = sceneGeometryFileName;
    tmp = strsplit(sceneGeometryIn,filesep);
    file = tmp{end};
    path = strcat(tmp(1:end-1));
end

% Load the selected sceneGeometry file
dataLoad=load(sceneGeometryIn);
sceneGeometrySource=dataLoad.sceneGeometry;
clear dataLoad

% Derive file stem path for the files associated with the sceneGeometry
pathParts = strsplit(path,filesep);
fileParts = strsplit(file,'_sceneGeometry.mat');
fileStem = fileParts{1};

% Define an alignment method. The choices are:
%    'gaze': If the sceneGeometry video was acquired while the subject was
%       performing a fixation task, and the acquisition includes an initial
%       period during which the subject fixated a centrally located target
%       on the screen, then align the pupil centers for images
%       corresponding to the eye oriented at a [0, 0] gaze position. This
%       is valid only for Session 2 TOME data collected after 01/30/2018.
%    'shape': Find in the acquisition video an image in which the pupil
%       aspect ratio and tilt is most similar to that seen while the
%       subject was at gaze position [0, 0].

% Use shape for all datasets except Session 2 data collected on or after
% 01/30/2018.

noFixSubjects = {'3001', '3002', '3003', '3004', '3005', '3007', '3008', '3009', '3011', '3013'};
if any(strcmp(pathParts,'session2_spatialStimuli')) && ~any(contains(pathParts,noFixSubjects))
    alignMethod = 'gaze';
    msg='Using the gaze alignment method';
    fprintf([msg '\n']);
    fprintf(fid, '%s: %s\n', datestr(now, 0), msg);
else
    alignMethod = 'shape';
    msg='Using the pupil shape alignment method';
    fprintf([msg '\n']);
    fprintf(fid, '%s: %s\n', datestr(now, 0), msg);
end

% Find the scene geometry video frames that correspond to fixation at the
% [0 0] screen position, and the eye at a neutral (zero azimuth and
% elevation) position
fileIn = fullfile(path,[fileStem,'_pupil.mat']);
load(fileIn,'pupilData');

% Obtain the eye rotation values from the pupilData, and convert these into
% gaze position on the screen.
eyeRotation = pupilData.radiusSmoothed.eyePoses.values(:,1:2)';
gazePosition = (sceneGeometrySource.screenPosition.R * eyeRotation + sceneGeometrySource.screenPosition.fixationAngles(1:2)')';

% Find the minimum fixation error threshold that results in a run of
% consecutive frames at fixation of the target length.
targetLength = 30;
runStarts = @(thresh) find(diff([0,(sqrt(sum(gazePosition.^2,2)) < thresh)',0]==1));
pullStartIndices = @(vec) vec(1:2:end-1);
pullRunLengths = @(vec) vec(2:2:end)-pullStartIndices(vec);
myObj = @(thresh) targetLength - max( pullRunLengths(runStarts(thresh)) );
threshVal = fzero(myObj,0.5);

% Find the start point of this run of frames
runLengths = pullRunLengths(runStarts(threshVal));
runIndices = pullStartIndices(runStarts(threshVal));
runLength = targetLength-myObj(threshVal);
startIndex = runIndices(runLengths == runLength);

% Load in the median image from the period of fixation for the
% sceneGeometry file. This is the "fixed" frame.
videoInFileName = fullfile(path,[fileStem '_gray.avi']);
fixedFrame = makeMedianVideoImage(videoInFileName,'startFrame',startIndex,'nFrames',runLength,'chunkSizeSecs',1/60);

% Find the median pupil center for these frames
fixFramePupilCenterFixation = [ ...
    nanmedian(pupilData.radiusSmoothed.ellipses.values(startIndex:startIndex+runLength,1)), ...
    nanmedian(pupilData.radiusSmoothed.ellipses.values(startIndex:startIndex+runLength,2)) ];

% Find the median theta and rho value for these frames (SEE:
% csaEllipseError)
fixFramePupilRhoShape = nanmedian(pupilData.radiusSmoothed.ellipses.values(startIndex:startIndex+runLength,4));
fixFramePupilRhoShape = 1-sqrt(1-fixFramePupilRhoShape^2);

fixFramePupilThetaShape = nanmedian(pupilData.radiusSmoothed.ellipses.values(startIndex:startIndex+runLength,5));
fixFramePupilThetaShape = fixFramePupilThetaShape*2;

% Load the perimter
perimeterFileName = fullfile(path,[fileStem '_correctedPerimeter.mat']);
load(perimeterFileName,'perimeter');

% Find the frame with the lowest ellipse RMSE during this period
rmseVals = pupilData.radiusSmoothed.ellipses.RMSE(startIndex:startIndex+runLength);
bestFrame = startIndex + find(rmseVals == min(rmseVals)) - 1;
Xpf = perimeter.data{bestFrame}.Xp;
Ypf = perimeter.data{bestFrame}.Yp;

% Get the camera offset point
cameraOffsetPoint = [sceneGeometrySource.cameraIntrinsic.matrix(1,3), ...
    sceneGeometrySource.cameraIntrinsic.matrix(2,3)];

% Find the pupil center for the eye model in the fixed image
eyePose = [0 0 0 3];
pupilEllipse = pupilProjection_fwd(eyePose,sceneGeometrySource);
fixFramePupilCenterZeroPos = pupilEllipse(1:2);

% Get a list of all gray.avi videos in this directory
fileList = dir(fullfile(path,'*_gray.avi'));

% Exclude the video that is the source of the fixed image
keep=cellfun(@(x) ~strcmp(x,[fileStem '_gray.avi']),extractfield(fileList,'name'));
fileList = fileList(keep);

% Ask the operator which of the videos we wish to adjust
fprintf('\n\nSelect the acquisition to adjust:\n')
for pp=1:length(fileList)
    optionName=['\t' num2str(pp) '. ' fileList(pp).name '\n'];
    fprintf(optionName);
end
fprintf('\nYou can enter a single acquisition number (e.g. 4),\n  a range defined with a colon (e.g. 4:7),\n  or a list within square brackets (e.g., [4 5 7]):\n')
choice = input('\nYour choice: ','s');
fileList = fileList(eval(choice));

% Create a figure
figHandle = figure();
imshow(fixedFrame,[],'Border','tight');
ax = gca;
ax.Toolbar = [];
hold on
text(20,30,'FIXED', 'Color', 'g','Fontsize',16);


% Provide some instructions for the operator
fprintf('Adjust horizontal and vertical camera translation with the arrow keys.\n');
fprintf('Adjust camera torsion with j and k.\n');
fprintf('Switch between moving and fixed image by pressing a.\n');
fprintf('Turn on and off perimeter display with p.\n');
fprintf('Turn on and off model display with m.\n');
fprintf('Press esc to exit.\n\n');
fprintf([path '\n']);

% Define a blank frame that we will need during display
blankFrame = ones(size(fixedFrame))*128;

% Loop over the selected acquisitions
for ff=1:length(fileList)
    
    % Load the timebase, pupilData, and perimeter for this acquisition
    acqFileStem = strsplit(fileList(ff).name,'_gray.avi');
    acqFileStem = acqFileStem{1};
    timebaseFileName = fullfile(path,[acqFileStem '_timebase.mat']);
    load(timebaseFileName,'timebase');
    pupilFileName =  fullfile(path,[acqFileStem '_pupil.mat']);
    load(pupilFileName,'pupilData');
    perimeterFileName =  fullfile(path,[acqFileStem '_correctedPerimeter.mat']);
    load(perimeterFileName,'perimeter');
    
    % Identify the startFrame, which is the time point at which the fMRI
    % acquisition began
    [~, startFrame] = min(abs(timebase.values));
    
    switch alignMethod
        case 'gaze'
            % Find the period of 30 frames prior to the start of the scan
            % when the eye was in the most consistent position, and closest
            % to the median position
            gazeX = pupilData.initial.ellipses.values(1:startFrame,1);
            gazeY = pupilData.initial.ellipses.values(1:startFrame,2);
            medianX = nanmedian(gazeX);
            medianY = nanmedian(gazeY);
            gazePosition = [gazeX-medianX; gazeY-medianY];
            
            targetLength = 30;
            runStarts = @(thresh) find(diff([0,(sqrt(sum(gazePosition.^2,2)) < thresh)',0]==1));
            pullStartIndices = @(vec) vec(1:2:end-1);
            pullRunLengths = @(vec) vec(2:2:end)-pullStartIndices(vec);
            myObj = @(thresh) targetLength - max( pullRunLengths(runStarts(thresh)) );
            threshVal = fzero(myObj,0.5);
        case 'shape'
            % Find the period of 30 frames after the start of the scan when
            % the pupil has a shape most similar to the shape from the
            % sceneGeometry file for gaze [0 0]
            rho = pupilData.initial.ellipses.values(startFrame:end,4);
            rho = 1-sqrt(1-rho.^2);
            theta = pupilData.initial.ellipses.values(startFrame:end,5);
            theta = theta.*2;
            
            shapeError = ...
                sqrt(fixFramePupilRhoShape^2 + rho.^2 - 2*fixFramePupilRhoShape.*rho.*cos(fixFramePupilThetaShape-theta))./2;
            
            targetLength = 30;
            runStarts = @(thresh) find(diff([0,(sqrt(sum(shapeError.^2,2)) < thresh)',0]==1));
            pullStartIndices = @(vec) vec(1:2:end-1);
            pullRunLengths = @(vec) vec(2:2:end)-pullStartIndices(vec);
            % This min([1e6 obj]) trick is to handle the objective otherwise
            % returning an empty value for a threshold of zero.
            myObj = @(thresh) min([1e6, targetLength - max( pullRunLengths(runStarts(thresh)) )]);
            threshVal = fzero(myObj,0.05);
            
    end
    
    % Find the start point of this run of frames
    runLengths = pullRunLengths(runStarts(threshVal));
    runIndices = pullStartIndices(runStarts(threshVal));
    runLength = targetLength-myObj(threshVal);
    startIndex = runIndices(runLengths == runLength);
    
    % Find the frame with the lowest ellipse RMSE during this period
    rmseVals = pupilData.initial.ellipses.RMSE(startIndex:startIndex+runLength);
    bestFrame = startIndex + find(rmseVals == min(rmseVals)) - 1;
    Xpm = perimeter.data{bestFrame}.Xp;
    Ypm = perimeter.data{bestFrame}.Yp;
    
    % If the fixationAngles in sceneGeometry is non-zero, then it may be
    % the case that this video to-be-adjusted has the property that the
    % subject was fixating the same screen location (0 0) at the start of
    % the scan. Make a guess at the x translation based upon this.
    if sceneGeometrySource.screenPosition.fixationAngles(1) ~= 0
        % Get the pupil center for the fames from the moving video
        movingFramePupilCenterFixation = [ ...
            nanmedian(pupilData.initial.ellipses.values(startIndex:startIndex+runLength,1)), ...
            nanmedian(pupilData.initial.ellipses.values(startIndex:startIndex+runLength,2)) ];
        % The adjustment is the difference in pupil centers from the fixed
        % and moving videos
        x = fixFramePupilCenterFixation - movingFramePupilCenterFixation;
    else
        movingFramePupilCenterFixation = [0 0];
        x = [0 0];
    end
    
    % Define the video file name
    videoInFileName = fullfile(path,fileList(ff).name);
    
    % Load the moving frame
    movingFrame = makeMedianVideoImage(videoInFileName,'startFrame',startIndex,'nFrames',runLength,'chunkSizeSecs',1/60);
    
    % Report which video we are working on
    msg=fileList(ff).name;
    fprintf(['\n' msg ':\n']);
    
    % Prepare for the loop
    showMoving = true;
    showPerimeter=false;
    showModel=false;
    annotHandle = [];
    torsion = 0;
    stillWorking = true;
    
    % Enter the while stillWorking loop
    while stillWorking
        
        % Prepare to update the image
        hold off
        
        if showMoving
            % Work with the moving frame
            displayImage = movingFrame;
            % Embed the movingFrame within a larger image that is padded
            % with mid-point background values
            padVals = round(size(movingFrame)./2);
            displayImagePad = zeros(size(movingFrame)+padVals.*2)+125;
            displayImagePad(padVals(1)+1:padVals(1)+size(movingFrame,1), ...
                padVals(2)+1:padVals(2)+size(movingFrame,2) ) = displayImage;
            displayImage = displayImagePad;
            % Apply the x and y translation
            displayImage = imtranslate(displayImage,x,'method','cubic');
            % Crop out the padding
            displayImage = displayImage(padVals(1)+1:padVals(1)+size(movingFrame,1), ...
                padVals(2)+1:padVals(2)+size(movingFrame,2));
            % Rotate the image
            displayImage = imrotateAround(displayImage, cameraOffsetPoint(2), cameraOffsetPoint(1), -torsion, 'bicubic');
            
            % Update the perimeter points
            [Xpa, Ypa] = updatePerimeter(Xpm,Ypm,x,-torsion,cameraOffsetPoint);
            
            % Display the perimeter points
            if showPerimeter
                idx = sub2ind(size(displayImage),round(Ypa),round(Xpa));
                displayImage(idx)=255;
            end
            
            % Display the image
            movingImHandle = imshow(displayImage,[],'Border','tight');
            ax = gca;
            ax.Toolbar = [];
            hold on
            text(20,30,'MOVING', 'Color', 'r','Fontsize',16);
            
        else
            displayImage = fixedFrame;
            
            % Update the perimeter
            Xpa = Xpf; Ypa = Ypf;
            
            % Display the perimeter points
            if showPerimeter
                idx = sub2ind(size(displayImage),round(Ypa),round(Xpa));
                displayImage(idx)=255;
            end
            
            fixedImHandle = imshow(displayImage,[],'Border','tight');
            ax = gca;
            ax.Toolbar = [];
            hold on
            text(20,30,'FIXED', 'Color', 'g','Fontsize',16);
        end
        
        % Show the eye model
        if showModel
            % Let the user know this will take a few seconds
            text_str = 'Updating model...';
            annotHandle = addAnnotation(text_str);
            % Obtain the eye pose from the adjusted perimeter
            eyePose = eyePoseEllipseFit(Xpa, Ypa, sceneGeometrySource);
            % Render the eye model
            renderEyePose(eyePose, sceneGeometrySource, ...
                'newFigure', false, 'visible', true, ...
                'showAzimuthPlane', true, ...
                'modelEyeLabelNames', {'retina' 'pupilEllipse' 'cornea'}, ...
                'modelEyePlotColors', {'.w' '-g' '.y'}, ...
                'modelEyeSymbolSizeScaler',1.5,...
                'modelEyeAlpha', 0.25);
            hold on
            % Remove the updating annotation
            delete(annotHandle);
        end
        
        % Add a marker for the camera CoP
        plot(cameraOffsetPoint(1),cameraOffsetPoint(2),'+c');
                
        keyAction = waitforbuttonpress;
        if keyAction
            keyChoiceValue = double(get(gcf,'CurrentCharacter'));
            switch keyChoiceValue
                case 28
                    x(1)=x(1)-1;
                case 29
                    x(1)=x(1)+1;
                case 30
                    x(2)=x(2)-1;
                case 31
                    x(2)=x(2)+1;
                case 97
                    showMoving = ~showMoving;
                case 106
                    torsion = torsion - 1;
                case 107
                    torsion = torsion + 1;
                case 112
                    showPerimeter = ~showPerimeter;
                case 109
                    showModel = ~showModel;
                case 27
                    text_str = 'finishing...';
                    annotHandle = addAnnotation(text_str);
                    stillWorking = false;
                otherwise
                    text_str = 'unrecognized command';
            end
            
        end
    end
    
    % Obtain the eye pose from the adjusted perimeter
    [Xpa, Ypa] = updatePerimeter(Xpm,Ypm,x,-torsion,cameraOffsetPoint);
    eyePose = eyePoseEllipseFit(Xpa, Ypa, sceneGeometrySource);

    % Calculate the updated camera rotation
    newCameraTorsion = sceneGeometrySource.cameraPosition.torsion - torsion;

    % Update the sceneGeometry torsion
    sceneGeometryAdjusted = sceneGeometrySource;
    sceneGeometryAdjusted.cameraPosition.torsion = newCameraTorsion;    

    % Find the change in the extrinsic camera translation needed to shift
    % the eye model the observed number of pixels
    p = calcCameraTranslationPixels(sceneGeometryAdjusted,eyePose,x);
        
    % Update the sceneGeometry translation
    sceneGeometryAdjusted.cameraPosition.translation = p;
    
    % Obtain the eye pose for the adjusted sceneGeometry
    eyePose = eyePoseEllipseFit(Xpm, Ypm, sceneGeometryAdjusted);
    
    % Report the values
    msg = sprintf('     adjustedCameraPositionTranslation [x; y; z] = [%2.3f; %2.3f; %2.3f]',p);
    fprintf(fid, '%s: %s\n', datestr(now, 0), msg);
    fprintf([msg '\n']);
    msg = sprintf('     adjustedCameraPositionTorsion [deg] = %2.3f',newCameraTorsion);
    fprintf(fid, '%s: %s\n', datestr(now, 0), msg);
    fprintf([msg '\n']);
    msg = sprintf('     adjustedFixationAngles [azi, ele, tor] = [%2.3f; %2.3f; %2.3f]',-eyePose(1:3));
    fprintf(fid, '%s: %s\n', datestr(now, 0), msg);
    fprintf([msg '\n']);
    
    % Show the original moving image and the adjusted sceneGeometry
    verifySolution = true;
    if verifySolution
        displayImage = movingFrame;
        idx = sub2ind(size(displayImage),round(Ypm),round(Xpm));
        displayImage(idx)=255;
        imshow(displayImage,[],'Border','tight');
        ax = gca;
        ax.Toolbar = [];
        hold on
        text(20,30,'MOVING - orig', 'Color', 'r','Fontsize',16);
        text_str = 'Updating model...';
        annotHandle = addAnnotation(text_str);
        renderEyePose(eyePose, sceneGeometryAdjusted, ...
            'newFigure', false, 'visible', true, ...
            'showAzimuthPlane', true, ...
            'modelEyeLabelNames', {'retina' 'pupilEllipse' 'cornea'}, ...
            'modelEyePlotColors', {'.w' '-g' '.y'}, ...
            'modelEyeSymbolSizeScaler',1.5,...
            'modelEyeAlpha', 0.25);
        delete(annotHandle);
        text_str = 'Adjusted geometry';
        annotHandle = addAnnotation(text_str);
        pause
    end
    
    % Revert to the fixed image
    delete(annotHandle)
    hold off
    fixedImHandle = imshow(fixedFrame,[],'Border','tight');
    ax = gca;
    ax.Toolbar = [];
    hold on
    text(20,30,'FIXED', 'Color', 'g','Fontsize',16);
    text_str = 'Loading...';
    annotHandle = addAnnotation(text_str);
    
end

close(figHandle);
fclose(fid);

end


%% LOCAL FUNCTIONS

function [Xp, Yp] = updatePerimeter(Xp,Yp,x,torsion,cameraOffsetPoint)

% Create a matrix of the perimeter points
v = [Xp';Yp'];

% Create the translation matrix
t = repmat([x(1); x(2)], 1, length(Xp));

% Translate the points
v = v+t;

% Set up the rotation matrix
center = repmat([cameraOffsetPoint(1); cameraOffsetPoint(2)], 1, length(Xp));
theta = deg2rad(-torsion);
R = [cos(theta) -sin(theta); sin(theta) cos(theta)];

% Apply the rotation
v = R*(v - center) + center;

% Extract the Xp and Yp vectors
Xp = v(1,:)';
Yp = v(2,:)';

end


function p = calcCameraTranslationPixels(sceneGeometrySource,eyePose,x)

% Find the change in the extrinsic camera translation needed to shift
% the eye model the observed number of pixels for an eye with zero rotation

p0 = sceneGeometrySource.cameraPosition.translation;
ub = sceneGeometrySource.cameraPosition.translation + [10; 10; 0];
lb = sceneGeometrySource.cameraPosition.translation - [10; 10; 0];
place = {'cameraPosition' 'translation'};
mySG = @(p) setfield(sceneGeometrySource,place{:},p);
pupilCenter = @(k) k(1:2);
targetPupilCenter = pupilCenter(pupilProjection_fwd(eyePose,sceneGeometrySource)) - x;
myError = @(p) norm(targetPupilCenter-pupilCenter(pupilProjection_fwd(eyePose,mySG(p))));
options = optimoptions(@fmincon,'Diagnostics','off','Display','off');
p = fmincon(myError,p0,[],[],[],[],lb,ub,[],options);
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

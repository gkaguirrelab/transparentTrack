% processVideosUnitTest
% 
% this script demonstrate the accuracy of the trasparentTrack engine on a
% "synthetic eye" video. The synthetic eye consists in a black ellipse and
% a white small dot (representing the pupil and the glint, respectively)
% moving on a gray background.
% 

% initialize the random number generator to the default seed
rng;

%% set path params
unitTestDir = '~/Desktop/unitTest';
if ~exist(unitTestDir,'dir')
    mkdir(unitTestDir)
end
syntheticEyeVideoName = fullfile(unitTestDir, 'syntheticEye.avi');
%% synthetic glint and pupil parameters

% 1. define synthetic data length
nFrames = 50;

% 2. define frame dimention in pixels
frameWidth = 320;
frameHeight = 240;

% 3. define glint position. As it is usually almost constant,  we pick a
% central location in the frame  and add some noise.
glint.X = random('Normal', frameWidth/2, 1, nFrames,1);
glint.Y = random('Normal', frameHeight/2, 1, nFrames,1);
glint.radius = 3;

% 4. define pupil behavior
% We'll have the pupil moving uniformely on a certain trajectory starting
% from the center of the frame, and insert a quick change of position
% halfway in the data to simulate a saccade.
pupil.X(1:nFrames/2) = linspace (frameWidth/2, frameWidth/2.5, nFrames/2 );
pupil.X(nFrames/2+1:nFrames) = linspace (frameWidth/2 + 10,  frameWidth/2.5+20,nFrames/2 );
pupil.X = pupil.X';
pupil.Y = linspace (frameHeight/2,frameHeight/2.2, nFrames )';

% For pupil size, assume a slowly decreasing pupil size over time.
pupil.width = linspace (82, 60, nFrames)';
pupil.height = linspace (70, 55, nFrames)';
pupil.theta = 0;


%%  generate and save video with these params

% create the video writer object
writerObj = VideoWriter(syntheticEyeVideoName, 'Uncompressed AVI');
writerObj.FrameRate = 60;

% open the video writer object
open(writerObj);

% initialize a gray empty frame
grayFrame=zeros(frameHeight, frameWidth, 3);
grayFrame(:,:,1) = .6;
grayFrame(:,:,2) =.6;
grayFrame(:,:,3) = .6;

% add pupil and glint to each frame
for ii =1:nFrames
    
    % get a certain  number of ellipse points that represent the pupil
    % boundary
    [xE,yE,~,~] = ellipse(3000, pupil.X(ii), pupil.Y(ii), pupil.width(ii)/2, pupil.height(ii)/2, pupil.theta);
    
    % build the "ellipse poligon" coordinates array
    poly = [];
    for jj = 1:length(xE)
        poly = [poly xE(jj) yE(jj)];
    end
    
    % add a filled black ellipse to the frame
    thisFrame = insertShape(grayFrame, 'filledPolygon', poly, 'Color', 'black','Opacity', 1);
    
    % add a circular white glint
    thisFrame = insertShape(thisFrame, 'filledCircle', [glint.X(ii) glint.Y(ii) glint.radius], 'Color', 'white','Opacity', 1);
    
    % store the frame
    outputVideo(:,:,:,ii)=thisFrame;
end

% Create a color map for the indexed video
cmap = [linspace(0,1,256)' linspace(0,1,256)' linspace(0,1,256)'];
cmap(1,:)=[1 0 0];
cmap(2,:)=[0 1 0];
cmap(3,:)=[0 0 1];
cmap(4,:)=[1 1 0];
cmap(5,:)=[0 1 1];
cmap(6,:)=[1 0 1];

% write the video
for ii=1:nFrames
    indexedFrame = rgb2ind(squeeze(outputVideo(:,:,:,ii)), cmap, 'nodither');
    writeVideo(writerObj,indexedFrame);
end


close (writerObj);


%% track synthetic eye
% define some file names
glintFileName = fullfile(unitTestDir, 'glint.mat');
perimeterFileName = fullfile(unitTestDir, 'perimeter.mat');
controlFileName = fullfile(unitTestDir, 'controlFile.csv');
correctedPerimeterFileName = fullfile(unitTestDir, 'correctedPerimeter.mat');
pupilFileName = fullfile(unitTestDir, 'pupil.mat');
finalFitVideoName = fullfile(unitTestDir, 'finalFit.avi');


% note: there is no need to use convertRawToGray in this case, because the
% video is digitally generated.

findGlint(syntheticEyeVideoName, glintFileName);

findPupilPerimeter(syntheticEyeVideoName, perimeterFileName);

makeControlFile(controlFileName, perimeterFileName, glintFileName, 'overwriteControlFile', true);

applyControlFile(perimeterFileName,controlFileName,correctedPerimeterFileName);

fitPupilPerimeter(correctedPerimeterFileName, pupilFileName);

makeFitVideo(syntheticEyeVideoName, finalFitVideoName, ...
    'glintFileName', glintFileName, 'perimeterFileName', correctedPerimeterFileName,...
    'pupilFileName', pupilFileName, 'whichFieldToPlot', 'pPosteriorMeanTransparent', ...
    'controlFileName',controlFileName);


%% compare tracking results with ground truth
% plot pupil
load(pupilFileName)
figure
subplot(2,1,1)
plot(pupilData.pPosteriorMeanTransparent(:,1))
hold on
plot(pupil.X)
xlabel('Frames')
ylabel('Pupil X [px]')
legend('fit', 'groundTruth')
title('X pupil position in pixels')
hold off

subplot(2,1,2)
plot(pupilData.pPosteriorMeanTransparent(:,2))
hold on
plot(pupil.Y)
xlabel('Frames')
ylabel('Pupil Y [px]')
legend('fit', 'groundTruth')
title('Y pupil position in pixels')
hold off

% plot glint
load(glintFileName)
figure
subplot(2,1,1)
plot(glintData.X)
hold on
plot(glint.X)
xlabel('Frames')
ylabel('Glint X [px]')
legend('fit', 'groundTruth')
title('X glint position in pixels')
hold off

subplot(2,1,2)
plot(glintData.Y)
hold on
plot(glint.Y)
xlabel('Frames')
ylabel('Glint Y [px]')
legend('fit', 'groundTruth')
title('Y glint position in pixels')
hold off


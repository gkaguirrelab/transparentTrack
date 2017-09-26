% transparentTrack unit test

% initialize the random number generator to the default seed
rng;

%% set path params
unitTestDir = '~/Desktop/unitTest';
if ~exist(unitTestDir,'dir')
    mkdir(unitTestDir)
end
syntheticEyeVideoName = fullfile(unitTestDir, 'syntheticEye.avi');
%% synthetic glint and pupil parameters

% define data length
nFrames = 50;
frameWidth = 320;
frameHeight = 240;
% the glint position is almost constant, so we just pick a value and add
% some noise
glint.X = random('Normal', frameWidth/2, 1, nFrames,1);
glint.Y = random('Normal', frameHeight/2, 1, nFrames,1);
glint.radius = 3;

% We have the pupil moving uniformely on a certain trajectory, and insert a
% quick change of position in the middle of the data
pupil.X(1:nFrames/2) = linspace (frameWidth/2, frameWidth/2.5, nFrames/2 );
pupil.X(nFrames/2+1:nFrames) = linspace (frameWidth/2 + 10,  frameWidth/2.5+20,nFrames/2 );
pupil.X = pupil.X';
pupil.Y = linspace (frameHeight/2,frameHeight/2.2, nFrames )';

% we assume a slowly decreasing pupil size over time
pupil.width = linspace (90, 40, nFrames)';
pupil.height = linspace (80, 35, nFrames)';
pupil.theta = 0;


%%  generate and save video with these params
% create the video writer with 1 fps
writerObj = VideoWriter(syntheticEyeVideoName, 'Uncompressed AVI');
writerObj.FrameRate = 60;

% open the video writer
open(writerObj);

grayFrame=zeros(frameHeight, frameWidth, 3);
grayFrame(:,:,1) = .6;
grayFrame(:,:,2) =.6;
grayFrame(:,:,3) = .6;
for ii =1:nFrames
    
    % add the pupil
    
    % get a certain  number of ellipse points
    [xE,yE,~,~] = ellipse(3000, pupil.X(ii), pupil.Y(ii), pupil.width(ii)/2, pupil.height(ii)/2, pupil.theta);
    
    % build the poligon coordinates
    poly = [];
    for jj = 1:length(xE)
        poly = [poly xE(jj) yE(jj)];
    end
    
    % build ellipse impicit equation
    thisFrame = insertShape(grayFrame, 'filledPolygon', poly, 'Color', 'black','Opacity', 1);
    
    
    % super impose the glint
    thisFrame = insertShape(thisFrame, 'filledCircle', [glint.X(ii) glint.Y(ii) glint.radius], 'Color', 'white','Opacity', 1);
    
    % store the frame
    outputVideo(:,:,:,ii)=thisFrame;
end

% write video

% Create a color map
cmap = [linspace(0,1,256)' linspace(0,1,256)' linspace(0,1,256)'];
cmap(1,:)=[1 0 0];
cmap(2,:)=[0 1 0];
cmap(3,:)=[0 0 1];
cmap(4,:)=[1 1 0];
cmap(5,:)=[0 1 1];
cmap(6,:)=[1 0 1];

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


% note: there is no need to use convertRawToGray

findGlint(syntheticEyeVideoName, glintFileName);

findPupilPerimeter(syntheticEyeVideoName, perimeterFileName);

makeControlFile(controlFileName, perimeterFileName, glintFileName);

applyControlFile(perimeterFileName,controlFileName,correctedPerimeterFileName);

fitPupilPerimeter(correctedPerimeterFileName, pupilFileName);

makeFitVideo(syntheticEyeVideoName, finalFitVideoName, ...
    'glintFileName', glintFileName, 'perimeterFileName', correctedPerimeterFileName,...
    'pupilFileName', pupilFileName, 'whichFieldToPlot', 'pPosteriorMeanTransparent', ...
    'controlFileName',controlFileName);


%% compare results
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


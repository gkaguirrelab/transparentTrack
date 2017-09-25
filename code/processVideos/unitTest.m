% transparentTrack unit test

%% set path params
systheticEyeVideoName = '~/Desktop/unitTest.avi';
%% synthetic glint and pupil parameters

% define data length
nFrames = 500;
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
 writerObj = VideoWriter(systheticEyeVideoName,'Grayscale AVI');
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
     
     thisFrame = rgb2gray(thisFrame);
     % write the frame
     writeVideo(writerObj, thisFrame);
 end
 
 close (writerObj);
 
 
 %% track synthetic eye
 
 
 
 %% compare results
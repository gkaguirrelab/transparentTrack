function frameArray = selectGazeCalFrames(pupilFileName, LTdatFileName, rawVidStartFileName, varargin)
% Select frames from gaze calibrtaion target fixation periods
%
% Syntax:
%  frameArray = selectGazeCalFrames(pupilFileName, LTdatFileName, rawVidStartFileName, varargin)
%
% Description:
%
% Inputs:
%	pupilFileName         - Full path to a pupilData file, a cell array
%                           of such paths, or a pupilData structure itself.
%                           If a single path, the pupilData file is loaded.
%                           If a cell array, the ellipse data from each
%                           pupilData file is loaded and concatenated.
%
% Optional key/value pairs (display and I/O):
%  'verbose'              - Logical. Default false.
%  'showPlots'            - Logical. Default false.
%
% Optional key/value pairs (analysis)
%  'fitLabel'             - Identifies the field in pupilData that contains
%                           the ellipse fit params for which the search
%                           will be conducted.
%
% Outputs
%	frameArray            - Vector.
%
% Examples:
%{
    pupilFileName = '/Users/aguirre/Dropbox (Aguirre-Brainard Lab)/TOME_processing/session2_spatialStimuli/TOME_3012/020317/EyeTracking/GazeCal02_pupil.mat';
    LTdatFileName = '/Users/aguirre/Dropbox (Aguirre-Brainard Lab)/TOME_data/session2_spatialStimuli/TOME_3012/020317/EyeTracking/GazeCal02_LTdat.mat';
    rawVidStartFileName = '/Users/aguirre/Dropbox (Aguirre-Brainard Lab)/TOME_data/session2_spatialStimuli/TOME_3012/020317/EyeTracking/GazeCal02_rawVidStart.mat';
    frameArray = selectGazeCalFrames(pupilFileName, LTdatFileName, rawVidStartFileName,'showPlots',true,'verbose',true);
%}


%% input parser
p = inputParser; p.KeepUnmatched = true;

% Required
p.addRequired('pupilFileName',@ischar);
p.addRequired('LTdatFileName',@ischar);
p.addRequired('rawVidStartFileName',@ischar);

% Optional display and I/O params
p.addParameter('verbose',false,@islogical);
p.addParameter('showPlots',false,@islogical);

% Optional analysis params
p.addParameter('fitLabel','initial',@ischar);

% parse
p.parse(pupilFileName, LTdatFileName, rawVidStartFileName, varargin{:})


%% Load pupil data
if ischar(pupilFileName)
    load(pupilFileName)
    ellipses = pupilData.(p.Results.fitLabel).ellipses.values;
    ellipseFitRMSE = pupilData.(p.Results.fitLabel).ellipses.RMSE;
end

%% Load live-track info files
LTGazeCalData=load(LTdatFileName);
dataLoad = load(rawVidStartFileName);
LTGazeCalData.rawVidStart = dataLoad.rawVidStart;
clear dataLoad

nTargets = size(LTGazeCalData.targets,1);

times = round((LTGazeCalData.dotTimes-LTGazeCalData.rawVidStart).*60);
xPos = pupilData.initial.ellipses.values(:,1);
yPos = pupilData.initial.ellipses.values(:,2);


% Define temporal support
support = 1:min([times(end)+5*60,size(xPos,1)]);
support=support(~isnan(xPos(support)));
support=support(~isnan(yPos(support)));

% Remove blinks
xVelocity = diff(xPos(support));
yVelocity = diff(yPos(support));
blinks = [0; (xVelocity>25) + (yVelocity>25)]';
for ii=-3:3
     blinks = blinks + circshift(blinks,ii);
end
blinks(blinks>0)=1;
support=support(~blinks);

% Normalize pupil position values
xPos = xPos-mean(xPos(support));
xPos = xPos./max(xPos(support));
yPos = yPos-mean(yPos(support));
yPos = yPos./max(yPos(support));

% Create the target vector
xTarget = zeros(size(xPos));
yTarget = zeros(size(yPos));
for ii=1:size(LTGazeCalData.targets)
    xTarget(times(ii):times(ii+1))=-sign(LTGazeCalData.targets(ii,1));
    yTarget(times(ii):times(ii+1))=sign(LTGazeCalData.targets(ii,2));
end


nonInteger = 1e6;
myObjX = @(x) 1-corr2(xPos(support),circshift(xTarget(support),round(x*nonInteger)));
myObjY = @(x) 1-corr2(yPos(support),circshift(yTarget(support),round(x*nonInteger)));
myObj = @(x) sqrt(mean([myObjX(x),myObjY(x)].^2));

[xShift,fVal]=fminsearch(myObj,0);
xShift = round(xShift*nonInteger);

% Find the lowest RMSE ellipse fit frame for each target
frames = times+xShift;
frameBoundEarly = frames(1:9)+round(diff(frames).*.35);
frameBoundLate = frames(1:9)+round(diff(frames).*.65);

for ii=1:nTargets
    supportStartIdx=find(frameBoundEarly(ii)<=support,1);
    supportEndIdx=find(frameBoundLate(ii)<=support,1);
    localSupport = support(supportStartIdx:supportEndIdx);
    [~,idx]=nanmin(ellipseFitRMSE(localSupport));
    frameArray(ii)=localSupport(idx);
end

% Plot the time series and the selected points
if p.Results.showPlots
figure
subplot(2,1,1)
plot(support/60,xPos(support),'-k');
hold on
plot(support/60,circshift(xTarget(support),xShift),'-b');
plot(frameArray/60,xPos(frameArray),'*r');
ylabel('xPos')
xlabel('time [sec]')
subplot(2,1,2)
plot(support/60,yPos(support),'-k');
hold on
plot(support/60,circshift(yTarget(support),xShift),'-b');
plot(frameArray/60,yPos(frameArray),'*r');
ylabel('yPos')
xlabel('time [sec]')

figure
plot(xPos(frameArray),yPos(frameArray),'bx');
end

if p.Results.verbose
    fprintf('Temporal offset of targets and pupil position: %0.0f frames (correlation: %0.2f) \n',xShift,1-fVal);
    outLine='ellipseArrayList: [ ';
    for ii=1:nTargets
        outLine = [outLine num2str(frameArray(ii))];
        if ii ~= nTargets
            outLine = [outLine ', '];
        end
    end
    outLine = [outLine ' ]\n'];
    fprintf(outLine);
end

end

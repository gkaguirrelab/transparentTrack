function [frameArray, fixationTargetArray] = selectGazeCalFrames(pupilFileName, LTdatFileName, rawVidStartFileName, varargin)
% Select frames from gaze calibration target fixation periods
%
% Syntax:
%  [frameArray, targetArray] = selectGazeCalFrames(pupilFileName, LTdatFileName, rawVidStartFileName, varargin)
%
% Description:
%   During gaze calibration acquisitions, the subject fixates a series of
%   targets of known position on a screen.
%
% Inputs:
%	pupilFileName         - Full path to a pupilData file.
%   LTdatFileName         - 
%   rawVidStartFileName   -
%
% Optional key/value pairs (display and I/O):
%  'verbose'              - Logical. Default false.
%  'showPlots'            - Logical. Default false.
%
% Optional key/value pairs (analysis)
%  'fitLabel'             - Identifies the field in pupilData that contains
%                           the ellipse fit params for which the search
%                           will be conducted.
%  'targetDeg'            - Scalar.
%
% Outputs
%	frameArray            - Vector.
%
% Examples:
%{
    pupilFileName = '~/Dropbox (Aguirre-Brainard Lab)/TOME_processing/session2_spatialStimuli/TOME_3012/020317/EyeTracking/GazeCal03_pupil.mat';
    LTdatFileName = '~/Dropbox (Aguirre-Brainard Lab)/TOME_data/session2_spatialStimuli/TOME_3012/020317/EyeTracking/GazeCal03_LTdat.mat';
    rawVidStartFileName = '~/Dropbox (Aguirre-Brainard Lab)/TOME_data/session2_spatialStimuli/TOME_3012/020317/EyeTracking/GazeCal03_rawVidStart.mat';
    [frameArray, fixationTargetArray] = selectGazeCalFrames(pupilFileName, LTdatFileName, rawVidStartFileName,'showPlots',true,'verbose',true);
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
p.addParameter('targetDeg',7,@isscalar);

% parse
p.parse(pupilFileName, LTdatFileName, rawVidStartFileName, varargin{:})


%% Load pupil data
if ischar(pupilFileName)
    load(pupilFileName)
    ellipses = pupilData.(p.Results.fitLabel).ellipses.values;
    ellipseFitRMSE = pupilData.(p.Results.fitLabel).ellipses.RMSE;
end

% Define temporal support
tmpDiff = diff(pupilData.timebase.values);
deltaT = tmpDiff(1);

%% Load live-track info files
LTGazeCalData=load(LTdatFileName);
dataLoad = load(rawVidStartFileName);
LTGazeCalData.rawVidStart = dataLoad.rawVidStart;
clear dataLoad

nTargets = size(LTGazeCalData.targets,1);

times = round((LTGazeCalData.dotTimes-LTGazeCalData.rawVidStart).*((1/deltaT)*1000));
xPos = pupilData.initial.ellipses.values(:,1);
yPos = pupilData.initial.ellipses.values(:,2);

support = 1:min([times(end)+5*(1/deltaT),size(xPos,1)]);
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
xTargetDegrees = -sign(LTGazeCalData.targets(:,1))*p.Results.targetDeg;
yTargetDegrees = -sign(LTGazeCalData.targets(:,2))*p.Results.targetDeg;
fixationTargetArray = [xTargetDegrees'; yTargetDegrees'];

% 
xTarget = zeros(size(xPos));
yTarget = zeros(size(yPos));
for ii=1:size(LTGazeCalData.targets)
    xTarget(times(ii):times(ii+1))= xTargetDegrees(ii);
    yTarget(times(ii):times(ii+1))= -yTargetDegrees(ii);
end


nonInteger = 1e6;
myObjX = @(x) 1-corr2(xPos(support),circshift(xTarget(support),round(x*nonInteger)));
myObjY = @(x) 1-corr2(yPos(support),circshift(yTarget(support),round(x*nonInteger)));
myObj = @(x) sqrt(mean([myObjX(x),myObjY(x)].^2));

[xShift,fVal]=fminsearch(myObj,0);
xShift = round(xShift*nonInteger);

% Find the lowest RMSE ellipse fit frame for each target
frames = times+xShift;
frameBoundEarly = frames(1:9)+round(diff(frames).*.15);
frameBoundLate = frames(1:9)+round(diff(frames).*.85);

for ii=1:nTargets
    supportStartIdx=find(frameBoundEarly(ii)<=support,1);
    supportEndIdx=find(frameBoundLate(ii)<=support,1);
    localSupport = support(supportStartIdx:supportEndIdx);    
    xPosMedian=weightedMedian(xPos(localSupport), ellipseFitRMSE(localSupport));
    yPosMedian=weightedMedian(yPos(localSupport), ellipseFitRMSE(localSupport));
    distanceFromMedian = ...
        sqrt(sum([(xPos(localSupport)-xPosMedian)';(yPos(localSupport)-yPosMedian)'].^2));
    [~,idx]=min(distanceFromMedian);
    frameArray(ii)=localSupport(idx);
end

% Plot the time series and the selected points
if p.Results.showPlots
figure
subplot(2,1,1)
plot(pupilData.timebase.values(support)./1000,xPos(support),'-k');
hold on
plot(pupilData.timebase.values(support)./1000,circshift(xTarget(support)./p.Results.targetDeg,xShift),'-b');
plot(pupilData.timebase.values(frameArray)./1000,xPos(frameArray),'*r');
ylabel('xPos')
xlabel('time [sec]')
subplot(2,1,2)
plot(pupilData.timebase.values(support)./1000,yPos(support),'-k');
hold on
plot(pupilData.timebase.values(support)./1000,circshift(yTarget(support)./p.Results.targetDeg,xShift),'-b');
plot(pupilData.timebase.values(frameArray)./1000,yPos(frameArray),'*r');
ylabel('yPos')
xlabel('time [sec]')

figure
plot(xPos(frameArray),yPos(frameArray),'bx');
end

if p.Results.verbose
    fprintf('Temporal offset of targets and pupil position: %0.0f frames (correlation: %0.2f) \n',xShift,1-fVal);
    outLine1='ellipseArrayList: [ ';
    outLine2='target array deg: [ ';
    outLine3=' ';
    for ii=1:nTargets
        outLine1 = [outLine1 num2str(frameArray(ii))];
        outLine2 = [outLine2 num2str(xTargetDegrees(ii))];
        outLine3 = [outLine3 num2str(yTargetDegrees(ii))];
        if ii ~= nTargets
            outLine1 = [outLine1 ', '];
            outLine2 = [outLine2 ', '];
            outLine3 = [outLine3 ', '];
        end
    end
    fprintf([outLine1 ' ]\n']);
    fprintf([outLine2 ' ;' outLine3 ']\n']);
end

end


function value = weightedMedian(data, weights)

weights = weights / sum(weights);

% sort the weights
[data,I] = sort(data);
weights = weights(I);

wd = weights.*data;
targetVal = sum(wd)/2;

value = [];      
j = 0;         
while isempty(value)
    j = j + 1;
    if targetVal >0
    if sum(wd(1:j)) >=targetVal
        value = data(j);    % value of the weighted median
    end
    else
    if sum(wd(1:j)) <=targetVal
        value = data(j);    % value of the weighted median
    end
    end
end

end
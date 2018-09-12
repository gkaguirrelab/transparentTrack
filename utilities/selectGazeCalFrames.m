function [frameArray, fixationTargetArray] = selectGazeCalFrames(pupilFileName, timebaseFileName, LTdatFileName, rawVidStartFileName, pupilCalInfoFileName, varargin)
% Select frames from gaze calibration target fixation periods
%
% Syntax:
%  [frameArray, targetArray] = selectGazeCalFrames(pupilFileName, LTdatFileName, rawVidStartFileName, varargin)
%
% Description:
%   During gaze calibration acquisitions, the subject fixates a series of
%   targets of known position on a screen. This routine uses information
%   regarding the timing and spatial position of fixation targets to
%   identify a set of video frames that best represent the position of the
%   eye during each of these fixation target periods.
%
% Inputs:
%	pupilFileName         - Full path to a pupilData file.
%	pupilFileName         - Full path to a timebase file.
%   LTdatFileName         -
%   rawVidStartFileName   - Full path to a rawVidStart file.
%   pupilCalInfoFileName  - Full path to a pupilCalInfo file. This file
%                           contains the information found in both the
%                           LTdat and the rawVidStart files. If defined,
%                           the values passed for LTdatFileName and
%                           rawVidStartFileName are ignored. A structure
%                           can also be passed for this file name, in which
%                           case the values in the structure are used
%                           instead of attempting to load the data from a
%                           file.
%
% Optional key/value pairs (display and I/O):
%  'verbose'              - Logical. Default false.
%  'showPlots'            - Logical. Default false.
%  'plotTitle'            - Char or empty.
%  'plotFileName'         - Char or empty. Full path that indicates where a
%                           diagnostic plot is to be saved. If empty, no
%                           file is saved.
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
    pupilFileName = '~/Dropbox (Aguirre-Brainard Lab)/TOME_processing/session2_spatialStimuli/TOME_3030/011018/EyeTracking/GazeCal01_pupil.mat';
    timebaseFileName = '~/Dropbox (Aguirre-Brainard Lab)/TOME_processing/session2_spatialStimuli/TOME_3030/011018/EyeTracking/GazeCal01_timebase.mat';
    gazeCalFramesDiagnosticPlot = '~/Dropbox (Aguirre-Brainard Lab)/TOME_processing/session2_spatialStimuli/TOME_3030/011018/EyeTracking/GazeCal01_fixFramesSelectPlot.pdf';
    LTdatFileName = '~/Dropbox (Aguirre-Brainard Lab)/TOME_data/session2_spatialStimuli/TOME_3030/011018/EyeTracking/GazeCal01_LTdat.mat';
    rawVidStartFileName = '~/Dropbox (Aguirre-Brainard Lab)/TOME_data/session2_spatialStimuli/TOME_3030/011018/EyeTracking/GazeCal01_rawVidStart.mat';
    pupilCalInfoFileName = '~/Dropbox (Aguirre-Brainard Lab)/TOME_data/session2_spatialStimuli/TOME_3030/011018/EyeTracking/GazeCal01_pupilCal_Info.mat';
    [frameArray, fixationTargetArray] = selectGazeCalFrames(pupilFileName, LTdatFileName, rawVidStartFileName,pupilCalInfoFileName,'plotFilename',gazeCalFramesDiagnosticPlot,'showPlot',true,'verbose',true);
%}


%% input parser
p = inputParser; p.KeepUnmatched = false;

% Required
p.addRequired('pupilFileName',@ischar);
p.addRequired('LTdatFileName',@ischar);
p.addRequired('rawVidStartFileName',@ischar);
p.addRequired('pupilCalInfoFileName',@(x)(ischar(x) || isstruct(x)));

% Optional display and I/O params
p.addParameter('verbose',false,@islogical);
p.addParameter('showPlot',false,@islogical);
p.addParameter('plotTitle','',@(x)(isempty(x) || ischar(x)));
p.addParameter('plotFileName','',@(x)(isempty(x) || ischar(x)));

% Optional analysis params
p.addParameter('fitLabel','initial',@ischar);
p.addParameter('targetDeg',7,@isscalar);

% parse
p.parse(pupilFileName, LTdatFileName, rawVidStartFileName, pupilCalInfoFileName,varargin{:})


%% Load pupil data
dataLoad = load(pupilFileName);
pupilData = dataLoad.pupilData;
clear dataLoad
ellipseFitRMSE = pupilData.(p.Results.fitLabel).ellipses.RMSE;

%% Load timebase
dataLoad = load(timebaseFileName);
timebase = dataLoad.timebase;
clear dataLoad

tmpDiff = diff(timebase.values);
deltaT = tmpDiff(1);

%% Load live-track info files
if isstruct(pupilCalInfoFileName)
    LTGazeCalData = pupilCalInfoFileName;
    xShift = 0;
else
    xShift = [];
    if exist(pupilCalInfoFileName, 'file') == 2
        dataLoad = load(pupilCalInfoFileName);
        LTGazeCalData.dotTimes = dataLoad.dotTimes;
        LTGazeCalData.targets = dataLoad.targets;
        LTGazeCalData.rawVidStart = dataLoad.rawVidStart.getSecsPre;
        clear dataLoad
    else
        if exist(LTdatFileName, 'file') == 2
            LTGazeCalData=load(LTdatFileName);
        else
            warning('There is no LTdata file for this acquisition; exiting')
            return
        end
        if exist(rawVidStartFileName, 'file') == 2
            dataLoad = load(rawVidStartFileName);
            LTGazeCalData.rawVidStart = dataLoad.rawVidStart;
            clear dataLoad
        else
            warning('There is no raw video start file for this acquisition; exiting')
            return
        end
    end
end

nTargets = size(LTGazeCalData.targets,1);

times = round((LTGazeCalData.dotTimes-LTGazeCalData.rawVidStart).*((1/deltaT)*1000));
xPosOriginal = pupilData.initial.ellipses.values(:,1);
yPosOriginal = pupilData.initial.ellipses.values(:,2);

support = 1:round(min([times(end)+5000*(1/deltaT),size(xPosOriginal,1)]));
support=support(~isnan(xPosOriginal(support)));
support=support(~isnan(yPosOriginal(support)));

% Remove blinks
xVelocity = diff(xPosOriginal(support));
yVelocity = diff(yPosOriginal(support));
blinks = [0; (xVelocity>25) + (yVelocity>25)]';
for ii=-3:3
    blinks = blinks + circshift(blinks,ii);
end
blinks(blinks>0)=1;
support=support(~blinks);

% Normalize pupil position values
xPos = xPosOriginal-mean(xPosOriginal(support));
xPos = xPos./max(xPos(support));
yPos = yPosOriginal-mean(yPosOriginal(support));
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

% Sometimes the target not tracked, and so was set as nans. Handles this
% here.
support=support(~isnan(xTarget(support)));

nonInteger = 1e6;
myObjX = @(x) 1-corr2(xPos(support),circshift(xTarget(support),round(x*nonInteger)));
myObjY = @(x) 1-corr2(yPos(support),circshift(yTarget(support),round(x*nonInteger)));
myObj = @(x) sqrt(mean([myObjX(x),myObjY(x)].^2));

if isempty(xShift)
    options = optimset('Display','off');
    [xShift,fVal]=fminsearch(myObj,0,options);
    if isnan(xShift)
        [xShift,fVal]=fminsearch(myObj,-20,options);
        if isnan(xShift)
            [xShift,fVal]=fminsearch(myObj,20,options);
        end
    end
    xShift = round(xShift*nonInteger);
else
    fVal = myObj(xShift);
end

% Find the lowest RMSE ellipse fit frame for each target
frames = times+xShift;
frameBoundEarly = frames(1:9)+round(diff(frames).*.15);
frameBoundLate = frames(1:9)+round(diff(frames).*.85);

for ii=1:nTargets
    if isnan(xTargetDegrees(ii))
        frameArray(ii)=nan;
    else
        supportStartIdx=find(frameBoundEarly(ii)<=support,1);
        supportEndIdx=find(frameBoundLate(ii)<=support,1);
        if isempty(supportEndIdx)
            if ii==9
                supportEndIdx=length(support);
            end
        end
        localSupport = support(supportStartIdx:supportEndIdx);
        xPosMedian=weightedMedian(xPos(localSupport), 1./ellipseFitRMSE(localSupport));
        yPosMedian=weightedMedian(yPos(localSupport), 1./ellipseFitRMSE(localSupport));
        distanceFromMedian = ...
            sqrt(sum([(xPos(localSupport)-xPosMedian)';(yPos(localSupport)-yPosMedian)'].^2));
        [~,idx]=min(distanceFromMedian);
        if isempty(idx)
            frameArray(ii)=nan;
        else
            frameArray(ii)=localSupport(idx);
        end
    end
end

% Remove any nan frames
notNanFrames = ~isnan(frameArray);
frameArray = frameArray(notNanFrames);
xTargetDegrees = xTargetDegrees(notNanFrames);
yTargetDegrees = yTargetDegrees(notNanFrames);
nTargets = length(xTargetDegrees);

if nTargets==0
    warning('No valid targets found; exiting')
    return    
end

% Plot the time series and the selected points
if p.Results.showPlot
    figHandle=figure('visible','on');
else
    figHandle=figure('visible','off');
end
% Configure the figure
set(gcf,'PaperOrientation','landscape');
set(figHandle, 'Units','inches')
height = 6;
width = 11;

% The last two parameters of 'Position' define the figure size
set(figHandle, 'Position',[25 5 width height],...
    'PaperSize',[width height],...
    'PaperPositionMode','auto',...
    'Color','w');

% Assemble some report text

outLineA = sprintf('Temporal offset: %0.0f frames; correlation: %0.2f',xShift,1-fVal);
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
outLineB = [outLine1 ' ]'];
outLineC = [outLine2 ' ;' outLine3 ']'];

if p.Results.verbose
    fprintf([outLineA ' \n']);
    fprintf([outLineB ' \n']);
    fprintf([outLineC ' \n']);
end


% Plot the data
if isempty(p.Results.plotTitle)
    nameParts = strsplit(pupilFileName,filesep);
    plotTitle = [nameParts{end-4} ' - ' nameParts{end-3} ' - ' nameParts{end-2} ' - ' nameParts{end}];
else
    plotTitle = p.Results.plotTitle;
end
subplot(3,4,[1 2])
plot(timebase.values(support)./1000,xPos(support),'-k');
hold on
plot(timebase.values(support)./1000,circshift(xTarget(support)./p.Results.targetDeg,xShift),'-b');
plot(timebase.values(frameArray)./1000,xPos(frameArray),'*r');
ylabel('xPos')
xlabel('time [sec]')
title(plotTitle,'Interpreter','none','HorizontalAlignment','left');

subplot(3,4,[5 6])
plot(timebase.values(support).*(1/deltaT),yPos(support),'-k');
hold on
plot(timebase.values(support).*(1/deltaT),circshift(yTarget(support)./p.Results.targetDeg,xShift),'-b');
plot(timebase.values(frameArray).*(1/deltaT),yPos(frameArray),'*r');
set(gca,'Ydir','reverse')
ylabel('yPos')
xlabel('time [frames]')

subplot(3,4,[3 4 7 8])
plot(xPosOriginal(frameArray),yPosOriginal(frameArray),'bx');
hold on
text(xPosOriginal(frameArray)+2,yPosOriginal(frameArray)+2,num2str(frameArray'));
ylim([min(yPosOriginal(frameArray)).*0.95 max(yPosOriginal(frameArray)).*1.05]);
xlim([min(xPosOriginal(frameArray)).*0.95 max(xPosOriginal(frameArray)).*1.05]);
set(gca,'Ydir','reverse')
axis equal

subplot(3,4,[9:12])
axis off
text(0.5,0.75,outLineA,'Units','normalized','HorizontalAlignment','center')
text(0.5,0.5,outLineB,'Units','normalized','HorizontalAlignment','center')
text(0.5,0.25,outLineC,'Units','normalized','HorizontalAlignment','center')


if ~isempty(p.Results.plotFileName)
    saveas(figHandle,p.Results.plotFileName)
end

if ~p.Results.showPlot
    close(figHandle)
end

end


function value = weightedMedian(data, weights)

weights = weights / sum(weights);
wd = weights.*data;

% sort the weighted data
[wd,I] = sort(wd);
data = data(I);

% Report the data value that corresponds to the median of the weighted
% values
value = data(find(wd>=median(wd),1));
end
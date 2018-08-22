function [] = selectDotTimesGUI(pupilFileName)

%% Load pupil data
if isempty(pupilFileName)
    [fileName, path] = uigetfile({'GazeCal*_pupil.mat'});
    if isempty(fileName)
        return
    end
    pupilFileName = [path, fileName];
end
load(pupilFileName)

xPosOriginal = pupilData.initial.ellipses.values(:,1);
yPosOriginal = pupilData.initial.ellipses.values(:,2);

% Define temporal support
if ~isfield(pupilData,'timebase')
    warning('This pupil data file does not include a timebase; exiting')
    return
end


% Define a suppot that cleans up the time series
support = 1:size(xPosOriginal,1);
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

xPosForPlotting = nan(size(xPosOriginal));
xPosForPlotting(support) = xPosOriginal(support);
yPosForPlotting = nan(size(yPosOriginal));
yPosForPlotting(support) = yPosOriginal(support);

figure
subplot(2,1,1)
xPosHandle = plot(xPosForPlotting,'-b');
hold on
subplot(2,1,2)
yPosHandle = plot(yPosForPlotting,'-b');
set(gca,'Ydir','reverse')
hold on

% Select border points
notDoneFlag = true;
frames = [];
while notDoneFlag
    [x1,y1] = ginput(1);
    if isempty(x1)
        notDoneFlag = false;
    else
    	plot(x1 ,y1, '+', 'Color', 'red');
        frames(end+1)=x1;
    end
end


nTargets = length(frames)-1;

frameBoundEarly = frames(1:nTargets);
frameBoundLate = frames(2:nTargets+1);

for ii=1:nTargets
        supportStartIdx=find(frameBoundEarly(ii)<=support,1);
        supportEndIdx=find(frameBoundLate(ii)<=support,1);
        localSupport = support(supportStartIdx:supportEndIdx);
        xPosMedian(ii)=median(xPosOriginal(localSupport));
        yPosMedian(ii)=median(yPosOriginal(localSupport));
end

if nTargets==9
    [~,targetRankX]=sort(xPosMedian);
    [~,targetRankY]=sort(yPosMedian);
    for ii=1:3
        idx = (1:3)+(ii-1)*3;
        xTarget(targetRankX(idx))=-(ii-2);
        yTarget(targetRankY(idx))=(ii-2);
    end    
else
    xClusters = kmeans(xPosMedian',3);
    for ii=1:3
        clusterMedian(ii)=median(xPosMedian(xClusters==ii));
    end
    [~,clusterRank]=sort(clusterMedian);
    for ii=1:3
        xTarget(xClusters==clusterRank(ii))=-(ii-2);
    end
    
    yClusters = kmeans(yPosMedian',3);
    for ii=1:3
        clusterMedian(ii)=median(yPosMedian(yClusters==ii));
    end
    [~,clusterRank]=sort(clusterMedian);
    for ii=1:3
        yTarget(yClusters==clusterRank(ii))=(ii-2);
    end
end

targets = [xTarget; yTarget]';

pupilCalInfo.dotTimes = pupilData.timebase.values(round(frames))./1000;
pupilCalInfo.targets = targets;
pupilCalInfo.rawVidStart = 0;

[a,b,c]=fileparts(pupilFileName);
plotFileName = fullfile(a,strrep([b c],'_pupil.mat','_fixFramesSelectPlot.pdf'));

selectGazeCalFrames(pupilFileName, '', '', pupilCalInfo ,'showPlot',true,'verbose',true,'plotFileName',plotFileName)

end

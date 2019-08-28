function selectDotTimesGUI(pupilFileName)
% A GUI interface to select the borders between gaze cal fixation trials
%
% Syntax:
%  selectDotTimesGUI(pupilFileName)
%
% Description:
%   During gaze calibration acquisitions, the subject fixates a series of
%   targets of known position on a screen. In some cases, the timing
%   of those fixation events is not known. This tool allows the user to
%   select a pupil.mat file that was derived from a gaze cal acquisition,
%   and then click on a plot to select the onset, transition times, and
%   final offset of fixation target events. The routine attempts to
%   determine the spatial positiom of those targets within a 3x3 grid, and
%   then passes the spatial position and timing of fixation events to the
%   routine selectGazeCalFrames.
%
% Inputs:
%	pupilFileName         - Full path to a pupilData file.
%
% Outputs
%	none
%
% Examples:
%{
    % ETTBSkip -- This is an interactive example.
    pupilFileName = '/Users/aguirre/Dropbox (Aguirre-Brainard Lab)/TOME_processing/session1_restAndStructure/TOME_3008/102116/EyeTracking/GazeCal_pupil.mat';
    selectDotTimesGUI(pupilFileName)
%}

%% Load pupil data
if isempty(pupilFileName)
    [fileName, path] = uigetfile({'GazeCal*_pupil.mat'});
    if isempty(fileName)
        return
    end
    pupilFileName = [path, fileName];
end

dataLoad=load(pupilFileName);
pupilData = dataLoad.pupilData;
clear dataLoad

% Extract the pupil X and Y positions from the initial ellipse fitting.
xPosOriginal = pupilData.initial.ellipses.values(:,1);
yPosOriginal = pupilData.initial.ellipses.values(:,2);

% Load the corresponding timebase file
[a,b,c]=fileparts(pupilFileName);
timebaseFileName = fullfile(a,[strrep(b,'_pupil','_timebase') c]);
dataLoad=load(timebaseFileName);
timebase = dataLoad.timebase;
clear dataLoad

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

% Plot the x and y pupil position time series
figure
subplot(2,1,1)
plot(xPosForPlotting,'-b');
ylabel('x pupil position [pixels]');
hold on
subplot(2,1,2)
plot(yPosForPlotting,'-b');
set(gca,'Ydir','reverse')
ylabel('y pupil position [pixels]');
xlabel('time [frames]');
hold on

% Select border points. Press return when done
notDoneFlag = true;
frames = [];
fprintf('Click on either the upper or lower plot to select the point of onset, \n');
fprintf('the transitions between, and the point of offset, of the fixation trials.\n');
fprintf('This will usually correspond to vertical lines in the pupil position plots.\n');
fprintf('To define the temporal boundaries of nine targets, place 10 markers.\n');
fprintf('Press return when finished placing markers.\n');
while notDoneFlag
    [x1,y1] = ginput(1);
    if isempty(x1)
        notDoneFlag = false;
    else
    	plot(x1 ,y1, '+', 'Color', 'red');
        frames(end+1)=x1;
    end
end

% Obtain the median vertical and horizontal pupil position for each of the
% periods marked by the defined frame bounds.
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

% Convert the medians into guesses regarding the target positions. This is
% well defined for the case of 9 targets (a 3x3 grid). If there are not 10
% boundary points (9 targets) the routine makes some guesses using
% clustering, but this is probably not trustworthy.
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

% Assemble the timing and spatial position information into a structure
pupilCalInfo.dotTimes = timebase.values(round(frames))./1000;
pupilCalInfo.targets = [xTarget; yTarget]';
pupilCalInfo.rawVidStart = 0;

% Define the filepath to be used to save a subsequent plot of gaze cal
% frames
[a,b,c]=fileparts(pupilFileName);
plotFileName = fullfile(a,strrep([b c],'_pupil.mat','_fixFramesSelectPlot.pdf'));

% Pass the assembled information to selectGazeCalFrames
selectGazeCalFrames(pupilFileName, '', '', pupilCalInfo ,'showPlot',true,'verbose',true,'plotFileName',plotFileName)

end

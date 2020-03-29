function saveModelFitPlot(obj,fileNameSuffix)


%% Obtain variables from the object
x = obj.x;
modelPupilEllipse = obj.modelPupilEllipse;
modelGlintCoord = obj.modelGlintCoord;
modelPoseGaze = obj.modelPoseGaze;
modelVecGaze = obj.modelVecGaze;
rawErrors = obj.rawErrors;


videoStemName = obj.videoStemName;
args = obj.args;
perimeter = args{1};
glintData = args{2};
ellipseRMSE = args{3};
gazeTargets = args{4};

plotFileName = [videoStemName '_sceneGeometry_modelFitPlot' fileNameSuffix '.pdf'];



% Prepare the figure
figHandle=figure('Visible','off');
set(gcf,'PaperOrientation','landscape');
set(figHandle, 'Units','inches')
height = 3;
width = 10;

% The last two parameters of 'Position' define the figure size
set(figHandle, 'Position',[25 5 width height],...
    'PaperSize',[width height],...
    'PaperPositionMode','auto',...
    'Color','w');


% We are going to have four sub-plots
nCols = 4;

% 4. Glint-pupil vec matching gaze targets
subplot(2,nCols,4)
plot(gazeTargets(1,:),gazeTargets(2,:),'ok'); hold on;
plot(modelVecGaze(1,:),modelVecGaze(2,:),'xr'); hold on;
ylim([-10 10])
axis equal
str = sprintf('Gaze vec [%2.2f]',rawErrors(4));
title(str);

% 3. EyePose matching gaze targets
subplot(2,nCols,3)
plot(gazeTargets(1,:),gazeTargets(2,:),'ok'); hold on;
plot(modelPoseGaze(1,:),modelPoseGaze(2,:),'xr'); hold on;
ylim([-10 10])
axis equal
str = sprintf('Gaze pose [%2.2f]',rawErrors(3));
title(str);

% 2. Glint fits
h1 = subplot(2,nCols,2);
plot(glintData.X,glintData.Y,'ok'); hold on;
plot(modelGlintCoord.X,modelGlintCoord.Y,'xr');
set(h1, 'Ydir', 'reverse')
axis equal
str = sprintf('Glint [%2.2f]',rawErrors(2));
title(str);

% 1. Perimeter fits
% Define a figure
hFig = figure( 'Visible', 'off');
dim = 150;
imshow(ones(dim,dim),'Border', 'tight');
drawnow;
hAxes = get(hFig,'CurrentAxes');
hold on;

% position the frames in the grid by their eye pose. First, figure out how
% many grid squares we have
nFrames = size(gazeTargets,2);
nGrid = ceil(sqrt(nFrames));

% Convert the gazeTargets array into an ordered list to store the frames
pos = gazeTargets ./ max(unique(abs(gazeTargets)));

% Reverse the sign of the elevation row, as positive eye movements
% correspond to the eye moving upwards in the image
pos(2,:) = -pos(2,:);

% Convet the positions into an ordered list
montageOrder = (pos(2,:)+nGrid-1).*nGrid+pos(1,:)-1;

% If we have a full supply of gazeTargets, use them to define the montage
% order. Otherwise, just go with the sortOrder
if all(~isnan(sum(gazeTargets)))
    montageOrder = montageOrder(sortOrder);
else
    montageOrder = 1:length(montageOrder);
end

% Show the frames
for ii = 1:length(ellipseRMSE)
    Xp = perimeter.data{ii}.Xp;
    meanXp = mean(Xp);
    Xp = Xp - meanXp + dim/2;
    Yp = perimeter.data{ii}.Yp;
    meanYp = mean(Yp);
    Yp = Yp - meanYp + dim/2;
    p1 = plot(hAxes,Xp,Yp,'.k');
    xlim([1 dim]);
    ylim([1 dim]);
    drawnow;
    hold on;
    pupilEllipseParams = modelPupilEllipse(ii,:);
    pupilEllipseParams(1) = pupilEllipseParams(1) - meanXp + dim/2;
    pupilEllipseParams(2) = pupilEllipseParams(2) - meanYp + dim/2;
    p2 = addTransparentEllipseToFigure(pupilEllipseParams,dim,dim,'red',1,hAxes);
    axis off;
    drawnow;
    thisFrame=getframe(hFig);
    framesToMontage(:,:,:,montageOrder(ii)) = thisFrame.cdata;
    delete(p1); delete(p2);
    drawnow;
end
close(hFig)
set(0, 'CurrentFigure', figHandle)
subplot(2,nCols,1)
montage(framesToMontage)
str = sprintf('Perimeter [%2.2f]',rawErrors(1));
title(str);

% Text label that indicates stage
xRange=get(gca,'XLim');
yRange=get(gca,'YLim');
drawnow


% Put an annotation for the x parameters at the
% bottom. Report params that hit a bound in red.
gcf;
axes('Position',[0 0 1 1],'Visible','off','Tag','subtitle');
dropboxBaseDir = getpref('eyeTrackTOMEAnalysis','dropboxBaseDir');
str = strrep(videoStemName,dropboxBaseDir,'');
ht = text(0.5,0.3,str,'Interpreter', 'none');
set(ht,'horizontalalignment','center','fontsize',10);
str = sprintf('Cornea curv joint, diff, tor, tilt, tip [$color-start$%2.2f$$color-end$$, $color-start$%2.2f$$color-end$$, $color-start$%2.2f$$color-end$$, $color-start$%2.2f$$color-end$$, $color-start$%2.2f$$color-end$$]; Rot center joint, diff [$color-start$%2.2f$$color-end$$, $color-start$%2.2f$$color-end$$]; primary pos [$color-start$%2.2f$$color-end$$, $color-start$%2.2f$$color-end$$]; Camera tor: $color-start$%2.1f$$color-end$$, trans: [$color-start$%2.1f$$color-end$$, $color-start$%2.1f$$color-end$$, $color-start$%2.1f$$color-end$$]',x);
tagIdx = strfind(str,'$color-start$');
fitAtBound = zeros(size(x));
for ii=1:length(fitAtBound)
    if fitAtBound(ii)
        str(tagIdx(ii):tagIdx(ii)+12) = '\color{red$$}';
    else
        str(tagIdx(ii):tagIdx(ii)+12) = '\color{black}';
    end
end
str = strrep(str,'$$color-end$$','\color{black}');
str = strrep(str,'$$','');
ht=text(.5,0.2,str);
set(ht,'horizontalalignment','center','fontsize',8);
str = sprintf('x = [ %2.3f, %2.3f, %2.3f, %2.3f, %2.3f, %2.3f, %2.3f, %2.3f, %2.3f, %2.3f, %2.3f, %2.3f, %2.3f ]',x);
ht=text(.5,0.1,str);set(ht,'horizontalalignment','center','fontsize',8);
drawnow

saveas(figHandle,plotFileName)

end

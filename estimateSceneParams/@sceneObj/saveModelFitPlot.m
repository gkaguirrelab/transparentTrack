function saveModelFitPlot(obj,fileNameSuffix)


%% Obtain variables from the object
x = obj.x;
model = obj.model;
modelPupilEllipse = obj.modelPupilEllipse;
modelGlintCoord = obj.modelGlintCoord;
modelPoseGaze = obj.modelPoseGaze;
modelVecGaze = obj.modelVecGaze;
rawErrors = obj.rawErrors;
perimeter = obj.perimeter;
glintDataX = obj.glintDataX;
glintDataY = obj.glintDataY;
frameSet = obj.frameSet;
gazeTargets = obj.gazeTargets;
videoStemName = obj.videoStemName;

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
plot(glintDataX,glintDataY,'ok'); hold on;
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

% If we don't have a full supply of gazeTargets, use the sortOrder to
% define the montage order.
if any(isnan(sum(gazeTargets)))
    [~, montageOrder] = sort(frameSet);
end

framesToMontage = (ones(dim,dim,3,nGrid^2));

% Show the frames
for ii = 1:length(frameSet)
    Xp = perimeter{ii}.Xp;
    meanXp = mean(Xp);
    Xp = Xp - meanXp + dim/2;
    Yp = perimeter{ii}.Yp;
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

% Annotate with the videoStemName
gcf;
axes('Position',[0 0 1 1],'Visible','off','Tag','subtitle');
dropboxBaseDir = getpref('eyeTrackTOMEAnalysis','dropboxBaseDir');
str = strrep(videoStemName,dropboxBaseDir,'');
ht = text(0.5,0.4,str,'Interpreter', 'none');
set(ht,'horizontalalignment','center','fontsize',10);

% Add the parameter values to the annotation
fields = {'head','eye','scene'};
for ii = 1:length(fields)
    str = [fields{ii} ' -- '];
    nParams = model.(fields{ii}).nParams;
    for pp = 1:nParams
        paramLabel = model.(fields{ii}).paramLabels{pp};
        str = [str paramLabel sprintf(': %2.2f',x(model.func.fieldParamIdx(fields{ii},paramLabel)))];
        if pp < nParams
            str = [str ', '];
        end
    end
    ht=text(.5,0.4 - 0.1*ii,str,'Interpreter', 'none');
    set(ht,'horizontalalignment','center','fontsize',10);
end
drawnow

saveas(figHandle,plotFileName)

end

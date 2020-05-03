function saveRelCameraPosPlot(obj,fileNameSuffix)
% Save a PDF of the relative camera position before and after adjustment
%
% Syntax:
%  obj.saveRelCameraPosPlot(fileNameSuffix)
%
% Description:
%   The search can adjust the relationship between the relative camera
%   position dimensions. This function saves a plot of the relative camera
%   position before and after adjustment
%
% Inputs:
%   fileNameSuffix        - Char vector. A suffix for the saved file name.
%                           The calling function might pass '1', '2', etc
%                           to label different stages of the evolution of
%                           the scene object.
%
% Outputs:
%   none
%


%% Exit if we do not have a relative camera position
if ~exist([obj.videoStemName '_relativeCameraPosition.mat'], 'file') == 2
    return
end


%% Setup variables

% Obtain variables from the object
x = obj.x;
model = obj.model;
videoStemName = obj.videoStemName;
origRelCamPos = obj.origRelCamPos;
relCamPos = obj.relCamPos;
nFrames = size(origRelCamPos,2);

% The name of the plot figure to be saved.
plotFileName = [videoStemName '_sceneGeometry_relCamPos' fileNameSuffix '.pdf'];


%% Prepare the figure

% Create a non-visible figure to hold the plot
figHandle=figure('Visible','off');
set(gcf,'PaperOrientation','landscape');
set(figHandle, 'Units','inches')
height = 6;
width = 10;

% The last two parameters of 'Position' define the figure size
set(figHandle, 'Position',[25 5 width height],...
    'PaperSize',[width height],...
    'PaperPositionMode','auto',...
    'Color','w');

subplot(2,1,1)
plot(origRelCamPos(1,:));
hold on
plot(origRelCamPos(2,:));
plot(origRelCamPos(3,:));
ylim([-5 5]);
xlim([0 ceil(nFrames/1000)*1000]);
legend({'+right','+up','+further'},'Location','eastoutside')
tLine1 = 'Relative camera position (initial)';
[~,tLine2] = fileparts(videoStemName);
tString = {tLine1,tLine2};
title(tString,'Interpreter','none');
xlabel('time [frames]');
ylabel('translation [mm]');
box off
legend boxoff

subplot(2,1,2)
plot(relCamPos(1,:));
hold on
plot(relCamPos(2,:));
plot(relCamPos(3,:));
ylim([-5 5]);
xlim([0 ceil(nFrames/1000)*1000]);
legend({'+right','+up','+further'},'Location','eastoutside')
tLine1 = 'Relative camera position (adjusted)';
% Add the parameter values to the title
tLine2 = '';
nParams = model.head.nParams;
for pp = 1:nParams
    paramLabel = model.head.paramLabels{pp};
    tLine2 = [tLine2 paramLabel sprintf(': %2.2f',x(model.func.fieldParamIdx('head',paramLabel)))];
    if pp < nParams
        tLine2 = [tLine2 ', '];
    end
end
tString = {tLine1,tLine2};
title(tString,'Interpreter','none');
xlabel('time [frames]');
ylabel('translation [mm]');
box off
legend boxoff
drawnow

% Save and close the figure
saveas(figHandle,plotFileName)
close(figHandle)

end
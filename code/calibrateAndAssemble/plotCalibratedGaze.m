function [gazePlot] = plotCalibratedGaze(gazeStruct,varargin)

% plotCalibratedGaze(gazeStruct)

% plots the calibrated gaze returning the figure handle for later figure
% adjustments and saving. The gaze can be plotted as a scatter plot or a
% timeseries in screen coordinates (mm on screen) or in polar coordinates.

%% Parse vargin for options passed here
p = inputParser; p.KeepUnmatched = true;

% Required
p.addRequired('gazeStruct',@isstruct);

% Optional analysis parameters
p.addParameter('whichCoordSystem','screen', @ischar); % alternative: polar
p.addParameter('plotType','scatter', @ischar); % alternative: timeseries
p.addParameter('customTitle','', @ischar); 
p.addParameter('screenWidth',698.5, @isnumeric);
p.addParameter('screenHeight',393.7, @isnumeric);
p.addParameter('screenRatio',[16 9], @isnumeric);
p.addParameter('calibratedUnits','mm', @ischar);


% Optional display and I/O parameters
p.addParameter('verbosity','none', @ischar);

% Environment parameters
p.addParameter('tbSnapshot',[],@(x)(isempty(x) | isstruct(x)));
p.addParameter('timestamp',char(datetime('now')),@ischar);
p.addParameter('username',char(java.lang.System.getProperty('user.name')),@ischar);
p.addParameter('hostname',char(java.net.InetAddress.getLocalHost.getHostName),@ischar);

% parse
p.parse(gazeStruct, varargin{:})



%% plot data
gazePlot = figure;
if strcmp(p.Results.whichCoordSystem,'screen')
    switch p.Results.plotType
        case 'scatter'
            %put a marker in the center of the screen
            plot(0,0,'+k')
            hold on
            % plot the scatter
            scatter(gazeStruct.X,gazeStruct.Y,6, 'filled', ...
                'MarkerEdgeColor', 'none', ...
                'MarkerFaceColor', [255 0 0]/255, ...
                'MarkerFaceAlpha', 0.01);
            xlim([-p.Results.screenWidth/2,p.Results.screenWidth/2])
            ylim([-p.Results.screenHeight/2,p.Results.screenHeight/2])
            xlabel('Screen width in mm from the center')
            ylabel('Screen height in mm from the center')
            set(gca,'Ydir','reverse')
            pbaspect([p.Results.screenRatio 1])
            
        case 'timeseries'
            % plot X coordinate of gaze
            subplot(2,1,1)
            plot(gazeStruct.X)
            xlabel('Frame')
            ylabel('X position of the gaze')
            % plot Y coordinate of gaze
            subplot(2,1,2)
            plot(gazeStruct.Y)
            xlabel('Frame')
            ylabel('Y position of the gaze')
            title(['Timeseries plot of gaze in screen coordinates (' p.Results.calibratedUnits ' on screen)'])
    end
    if isempty(p.Results.customTitle)
        title([p.Results.plotType 'plot of gaze in screen coordinates (' p.Results.calibratedUnits ' on screen)'])
    else
        title(p.Results.customTitle)
    end
end

if strcmp(p.Results.whichCoordSystem,'polar')
    switch p.Results.plotType
        case 'scatter'
            polarscatter(deg2rad(gazeStruct.pol),gazeStruct.ecc,6, 'filled', ...
                'MarkerEdgeColor', 'none', ...
                'MarkerFaceColor', [255 0 0]/255, ...
                'MarkerFaceAlpha', 0.01);
            pax = gca;
            pax.ThetaDir = 'clockwise';
            pax.ThetaZeroLocation = 'top';
            pax.RLim = [0 10];
        case 'timeseries'
            % plot eccentricity of gaze
            subplot(2,1,1)
            plot(gazeStruct.ecc)
            xlabel('Frame')
            ylabel('Eccentricity')
            % plot polar angle of gaze
            subplot(2,1,2)
            plot(gazeStruct.pol)
            xlabel('Frame')
            ylabel('Polar Angle')
    end
    if isempty(p.Results.customTitle)
        title([p.Results.plotType 'plot of gaze in polar coordinates'])
    else
        title(p.Results.customTitle)
    end
end
end
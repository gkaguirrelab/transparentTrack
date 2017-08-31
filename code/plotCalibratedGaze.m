function plotCalibratedGaze(gazeStruct,varargin)

% plotCalibratedGaze(gazeStruct)

% got to write this

%% Parse vargin for options passed here
p = inputParser; p.KeepUnmatched = true;

% Required
p.addRequired('gazeStruct',@isstruct);

% Optional analysis parameters
p.addParameter('whichCoordSystem','both', @ischar); % alternative: polar, both
p.addParameter('plotType','scatter', @ischar); % alternative: timeseries
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
if strcmp(p.Results.whichCoordSystem,'both') || strcmp(p.Results.whichCoordSystem,'screen')
    switch p.Results.plotType
        case 'scatter'
            figure
            %put a marker in the center
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
            title(['Scatter plot of gaze in screen coordinates (' p.Results.calibratedUnits ' on screen)'])
        case 'timeseries'
            % got to write this
    end
end

if strcmp(p.Results.whichCoordSystem,'both') || strcmp(p.Results.whichCoordSystem,'polar')
    switch p.Results.plotType
        case 'scatter'
            figure
            polarscatter(deg2rad(gazeStruct.pol),gazeStruct.ecc,6, 'filled', ...
                'MarkerEdgeColor', 'none', ...
                'MarkerFaceColor', [255 0 0]/255, ...
                'MarkerFaceAlpha', 0.01);
            pax = gca;
            pax.ThetaDir = 'clockwise';
            pax.ThetaZeroLocation = 'top';
            pax.RLim = [0 10];
            title('Scatter plot of gaze in polar coordinates')
        case 'timeseries'
            % got to write this
    end
end
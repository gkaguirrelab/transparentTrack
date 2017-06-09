function [pupilSize,gaze, pupilError,pupilCut] = calcPupilGaze(params)

% Calculates the pupil size and gaze position from eye tracking data
%
%   Usage:
%       [pupilSize,gaze] = calcPupilGaze(params)
%
%   If params.trackType == 'LiveTrack' <default>
%
%       Required inputs:
%           params.scaleCalFile     - full path to scale calibration .mat file
%           params.gazeCalFile      - full path to scale calibration .mat file
%           params.eyeTrackFile     - full path to eye tracking .mat file
%
%   If params.trackType == 'trackPupil'
%
%       Required inputs:
%           params.scaleCalVideo    - full path to scale calibration .avi video file
%           params.calTargetFile    - full path to gaze calibration LTdat.mat file
%           params.gazeCalVideo     - full path to gaze calibration .mov/.avi video file
%           params.eyeTrackVideo    - full path to eye tracking .avi video file
%
%       Optional inputs:
%           params.scaleCalOutVideo - full path to output scale calibration video
%           params.scaleCalOutMat   - full path to output scale calibration .mat file
%           params.gazeCalOutVideo  - full path to output gaze calibration video
%           params.gazeCalOutMat    - full path to output gaze calibration .mat file
%           params.eyeTrackOutVideo - full path to output eye tracking video
%           params.eyeTrackOutMat   - full path to output eye tracking .mat file
%
%   Defaults:
%       params.scaleSize        - 5;        % calibration dot size (mm)
%       params.vidBuffer        - 0.25;     % proportion of the gazeCalVideo to crop
%       params.viewDist         - 1065;     % distance from eyes to screen (mm)
%       params.trackType        - 'LiveTrack' (other option is 'trackPupil')
%
%   Outputs:
%       pupilSize               - vector of pupil sizes (mm)
%       gaze.X                  - vector of gaze X coordinates (mm)
%       gaze.Y                  - vector of gaze Y coordinates (mm)
%       gaze.ecc                - vector of gaze eccentricity values (degrees visual angle)
%       gaze.pol                - vector of gaze polar angle values (degrees)
%
%   Note about blinks:
%
%       If params.trackType == 'LiveTrack'
%           User must filter these out post-hoc, as the LiveTrack assigns a
%           value to the pupil and glint no matter what
%
%       If params.trackType == 'trackPupil'
%           nans in the output data indicate a blink, or really any event
%           where both the pupil and glint could not be simultaneously
%           tracked
%
%       If params.trackType == 'Hybrid'
%           The input is a "rescaled" version of tracked data obtained with
%           track pupil, that will be processed as if they were LiveTrack
%           data.
%
%   Written by Andrew S Bock Oct 2016

%% set defaults
if ~isfield(params,'scaleSize')
    params.scaleSize = 5;
end
if ~isfield(params,'vidBuffer')
    params.vidBuffer = 0.25;
end
if ~isfield(params,'viewDist')
    params.viewDist = 1065;
end
if ~isfield(params,'trackType')
    params.trackType = 'LiveTrack';
end
%% Get the mm / pixel from calibration stick
switch params.trackType
    case 'trackPupil'
        scaleParams.inVideo     = params.scaleCalVideo;
        if isfield(params,'scaleCalOutVideo');
            scaleParams.outVideo  = params.scaleCalOutVideo;
        end
        if isfield(params,'scaleCalOutMat');
            scaleParams.outMat  = params.scaleCalOutMat;
        end
        scaleParams.pupilRange  = params.scaleSize * [3 10];
        scaleParams.threshVals  = [0.05 0.999]; % bin for pupil and glint, respectively
        scaleParams.pupilOnly   = 1;
        [scalePupil]            = trackPupil(scaleParams);
        mmPerPixel              = params.scaleSize / median(scalePupil.size);
    case 'LiveTrack'
        scaleCal                = load(params.scaleCalFile);
        [maxVal,maxInd]         = max(scaleCal.ScaleCal.pupilDiameterMmGroundTruth);
        mmPerPixel              = maxVal / median([scaleCal.ScaleCal.ReportRaw{maxInd}.PupilWidth_Ch01]);
        pxPerMm = scaleCal.ScaleCal.cameraUnitsToMmWidthMean;
    case 'Hybrid'
        if isfield(params, 'sizeConversionFactor')
            pxPerMm = params.sizeConversionFactor.sizeConversionFactor;
        else % just use livetrack data
            scaleCal                = load(params.scaleCalFile);
            [maxVal,maxInd]         = max(scaleCal.ScaleCal.pupilDiameterMmGroundTruth);
            mmPerPixel              = maxVal / median([scaleCal.ScaleCal.ReportRaw{maxInd}.PupilWidth_Ch01]);
            pxPerMm = scaleCal.ScaleCal.cameraUnitsToMmWidthMean;
        end
        
end
%% Get the transformation matrix for gaze
switch params.trackType
    case 'trackPupil'
        % Track the pupil and glint in the gaze calibration video
        gazeParams.inVideo      = params.gazeCalVideo;
        if isfield(params,'gazeCalOutVideo')
            gazeParams.outVideo = params.gazeCalOutVideo;
        end
        if isfield(params,'gazeCalOutMat')
            gazeParams.outMat = params.gazeCalOutMat;
        end
        [gazePupil,gazeGlint]   = trackPupil(gazeParams);
        % The video is manually stopped, so there is too much video at the end
        % after the dots are gone
        pad                     = round(params.vidBuffer*length(gazePupil.X));
        gazePupil.X             = gazePupil.X(1:end-pad);
        gazePupil.Y             = gazePupil.Y(1:end-pad);
        gazePupil.size          = gazePupil.size(1:end-pad);
        gazeGlint.X             = gazeGlint.X(1:end-pad);
        gazeGlint.Y             = gazeGlint.Y(1:end-pad);
        gazeGlint.size          = gazeGlint.size(1:end-pad);
        % Pull out pupil and glint values for each point
        nPoints                 = 20; % last nPoints
        prevPoints              = zeros(nPoints,2);
        prevDist                = nan(size(gazePupil.X));
        for i = 1:length(gazePupil.X)
            clear tmpDist
            for j = 1:nPoints
                tmpDist(j)      =  sqrt( (gazePupil.X(i) - prevPoints(j,1))^2 + ...
                    (gazePupil.Y(i) - prevPoints(j,2))^2 );
            end
            prevDist(i)         = nanmean(tmpDist);
            prevPoints          = [prevPoints(2:end,:);[gazePupil.X(i),gazePupil.Y(i)]];
        end
        % Thresh the distances
        % Look at this plot, confirm the clusters are accurate
        thresh                  = nPoints/10;
        k                       = 9;
        goodInd                 = prevDist<thresh;
        gazePupilX              = gazePupil.X(goodInd);
        gazePupilY              = gazePupil.Y(goodInd);
        gazeGlintX              = gazeGlint.X(goodInd);
        gazeGlintY              = gazeGlint.Y(goodInd);
        % Provide the starting points for the dot search
        Xs                      = [min(gazePupilX), min(gazePupilX) + ...
            (max(gazePupilX) - min(gazePupilX))/2, max(gazePupilX)];
        Ys = [min(gazePupilY), min(gazePupilY) + ...
            (max(gazePupilY) - min(gazePupilY))/2, max(gazePupilY)];
        dotMatrix = [...
            Xs(1),Ys(1); ...
            Xs(1),Ys(2); ...
            Xs(1),Ys(3); ...
            Xs(2),Ys(1); ...
            Xs(2),Ys(2); ...
            Xs(2),Ys(3); ...
            Xs(3),Ys(1); ...
            Xs(3),Ys(2); ...
            Xs(3),Ys(3); ...
            ];
        % cluster the data
        idx = kmeans([gazePupilX,gazePupilY],k,'Start',dotMatrix);
        % Plot the means
        gazePupilMeans          = nan(k,2);
        gazeGlintMeans          = nan(k,2);
        fullFigure;
        for i = 1:k
            subplot(3,3,i);
            plot(gazePupilX,gazePupilY,'.','MarkerSize',20);
            axis square;
            hold on;
            title(['Cluster = ' num2str(i)],'FontSize',20);
            gazePupilMeans(i,1) = mean(gazePupilX(idx==i));
            gazePupilMeans(i,2) = mean(gazePupilY(idx==i));
            gazeGlintMeans(i,1) = mean(gazeGlintX(idx==i));
            gazeGlintMeans(i,2) = mean(gazeGlintY(idx==i));
            plot(gazePupilMeans(i,1),gazePupilMeans(i,2),'.r','MarkerSize',10);
            plot(gazeGlintMeans(i,1),gazeGlintMeans(i,2),'.g','MarkerSize',10);
        end
        % Pull out the pupil, glint, and target values
        gazeCal                 = load(params.calTargetFile);
        targets                 = gazeCal.targets;
        % Re-order pupil and glint to match order of targets (only if using 'trackPupil' and above code)
        uTX = unique(targets(:,1));
        uTY = unique(targets(:,2));
        ct = 0;
        for i = 1:length(uTX)
            for j = 1:length(uTY)
                ct = ct + 1;
                targInd = find(targets(:,1) == uTX(i) & targets(:,2) == uTY(j));
                calParams.pupil.X(targInd) = gazePupilMeans(ct,1);
                calParams.pupil.Y(targInd) = gazePupilMeans(ct,2);
                calParams.glint.X(targInd) = gazeGlintMeans(ct,1);
                calParams.glint.Y(targInd) = gazeGlintMeans(ct,2);
            end
        end
        % Calculate the 'CalMat'
        calParams.targets.X     = targets(:,1); % mm on screen, screen center = 0
        calParams.targets.Y     = targets(:,2); % mm on screen, screen center = 0
        calParams.viewDist      = params.viewDist; % mm from screen
        % Calculate the adjustment factor
        calParams.rpc           = calcRpc(calParams);
        % Calculate the transformation matrix
        [calParams.calMat]      = calcCalMat(calParams);
        calGaze                 = calcGaze(calParams);
        figure;
        hold on;
        % plot each true and tracked target position. Red cross means target
        % position and blue means tracked gaze position.
        for i = 1:length(calGaze.X)
            plot(targets(i,1), targets(i,2),'rx');
            plot(calGaze.X(i),calGaze.Y(i),'bx');
            plot([targets(i,1) calGaze.X(i)], [targets(i,2) calGaze.Y(i)],'g');
        end
    case {'LiveTrack'}
        gazeCal                 = load(params.gazeCalFile);
        calParams.rpc           = gazeCal.Rpc;
        calParams.calMat        = gazeCal.CalMat;
    case {'Hybrid'}
        gazeCal                 = load(params.gazeCalFile);
        if isfield(gazeCal, 'calParams')
            calParams.rpc           = gazeCal.calParams.rpc;
            calParams.calMat        = gazeCal.calParams.calMat;
        else
            calParams.rpc           = gazeCal.Rpc;
            calParams.calMat        = gazeCal.CalMat;
        end
end
%% Get the pupil size and gaze location from an eye tracking video
switch params.trackType
    case 'trackPupil'
        eyeParams.inVideo       = params.eyeTrackVideo;
        if isfield(params,'eyeTrackOutVideo');
            eyeParams.outVideo  = params.eyeTrackOutVideo;
        end
        if isfield(params,'eyeTrackOutMat');
            eyeParams.outMat    = params.eyeTrackOutMat;
        end
        [eyePupil,eyeGlint]     = trackPupil(eyeParams);
        eyeParams.pupil.X       = eyePupil.X;
        eyeParams.pupil.Y       = eyePupil.Y;
        eyeParams.glint.X       = eyeGlint.X;
        eyeParams.glint.Y       = eyeGlint.Y;
        eyeParams.viewDist      = params.viewDist;
        eyeParams.rpc           = calParams.rpc;
        eyeParams.calMat        = calParams.calMat;
    case 'LiveTrack'
        eyeMat                 = load(params.eyeTrackFile);
        eyePupil.size           = [];
        eyeParams.pupil.X       = [];
        eyeParams.pupil.Y       = [];
        eyeParams.glint.X       = [];
        eyeParams.glint.Y       = [];
        for i = 1:length(eyeMat.Report)
            % pupil size
            eyePupil.size       = [eyePupil.size;...
                eyeMat.Report(i).PupilWidth_Ch01;eyeMat.Report(i).PupilWidth_Ch02];
            % pupil X
            eyeParams.pupil.X       = [eyeParams.pupil.X;...
                eyeMat.Report(i).PupilCameraX_Ch01;eyeMat.Report(i).PupilCameraX_Ch02];
            % pupil Y
            eyeParams.pupil.Y       = [eyeParams.pupil.Y;...
                eyeMat.Report(i).PupilCameraY_Ch01;eyeMat.Report(i).PupilCameraY_Ch02];
            % glint X
            eyeParams.glint.X       = [eyeParams.glint.X;...
                eyeMat.Report(i).Glint1CameraX_Ch01;eyeMat.Report(i).Glint1CameraX_Ch02];
            % glint Y
            eyeParams.glint.Y       = [eyeParams.glint.Y;...
                eyeMat.Report(i).Glint1CameraY_Ch01;eyeMat.Report(i).Glint1CameraY_Ch02];
        end
        eyeParams.viewDist      = params.viewDist;
        eyeParams.rpc           = calParams.rpc;
        eyeParams.calMat        = calParams.calMat;
    case 'Hybrid'
        eyeMat                  = load(params.eyeTrackFile);
        eyePupil.size           = eyeMat.pupil.size;
        eyeParams.pupil.X       = eyeMat.pupil.X;
        eyeParams.pupil.Y       = eyeMat.pupil.Y;
        eyeParams.glint.X       = eyeMat.glint.X;
        eyeParams.glint.Y       = eyeMat.glint.Y;
        eyeParams.viewDist      = params.viewDist;
        eyeParams.rpc           = calParams.rpc;
        eyeParams.calMat        = calParams.calMat;
end
eyeGaze                         = calcGaze(eyeParams);
%% Set outputs

% pupilSize                       = eyePupil.size * mmPerPixel;
pupilSize                       = eyePupil.size ./ pxPerMm;
gaze                            = eyeGaze;

% also, save out the errorMetric for the pupil
pupilError = nan(length(eyeMat.pupil.X),1);
pupilCut = nan(length(eyeMat.pupil.X),1);
for ii = 1:length(eyeMat.pupil.X)
    if eyeMat.pupil.flags.cutPupil(ii) == 1
        pupilError(ii) = eyeMat.pupil.cutDistanceErrorMetric(ii);
        pupilCut(ii) = 1;
    else
        pupilError(ii) = eyeMat.pupil.distanceErrorMetric(ii);
        pupilCut(ii) = 0;
    end
end
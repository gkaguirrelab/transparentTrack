% DEMO_eyeTracking
%
% Demonstrate the entire eyetracking analysis pipeline.
%
% A local sandbox folder named 'eyeTrackingDEMO' will be created on the
% desktop to replicate the dropbox environment of the real routine. Files
% will be downloaded from figshare and placed in the sandbox (about 7 GB).
%
% Make sure your machine is configured to work with ToolboxToolbox.
%
% Run-time on an average computer is about 5 minutes for 500 frames.
% Set nFrames to 'Inf' to process the entire video, which will take longer.
%
% Usage examples
% ==============
%
% DEMO_eyeTracking
%


%% hard coded parameters
nFrames = Inf; % number of frames to process (set to Inf to do all)
verbose = true; % Set to none to make the demo silent
TbTbToolboxName = 'transparentTrack';


%% set paths and make directories
% create test sandbox on desktop
sandboxDir = '~/Desktop/eyeTrackingDEMO';
if ~exist(sandboxDir,'dir')
    mkdir(sandboxDir)
end

% define path parameters
pathParams.dataSourceDirRoot = fullfile(sandboxDir,'TOME_data');
pathParams.dataOutputDirRoot = fullfile(sandboxDir,'TOME_processing');
pathParams.projectSubfolder = 'session2_spatialStimuli';
pathParams.eyeTrackingDir = 'EyeTracking';
pathParams.subjectID = 'TOME_3020';
pathParams.sessionDate = '050517';
pathParams.runName = 'GazeCal01';


%% TbTb configuration
% We will suppress the verbose output, but detect if there are deploy
% errors and if so stop execution
tbConfigResult=tbUse(TbTbToolboxName,'reset','full','verbose',false);
if sum(cellfun(@sum,extractfield(tbConfigResult, 'isOk')))~=length(tbConfigResult)
    error('There was a tb deploy error. Check the contents of tbConfigResult');
end
% We save a deployment snapshot. This variable is passed to the analysis
% pipeline and then saved with every output file, thereby documenting the
% system and software configuration at the time of execution.
tbSnapshot=tbDeploymentSnapshot(tbConfigResult,[],'verbose',false);
clear tbConfigResult


%% Prepare paths and directories
% define full paths for input and output
pathParams.dataSourceDirFull = fullfile(pathParams.dataSourceDirRoot, pathParams.projectSubfolder, ...
    pathParams.subjectID, pathParams.sessionDate, pathParams.eyeTrackingDir);
pathParams.dataOutputDirFull = fullfile(pathParams.dataOutputDirRoot, pathParams.projectSubfolder, ...
    pathParams.subjectID, pathParams.sessionDate, pathParams.eyeTrackingDir);

% Download the data if it is not already there
demoPackage = fullfile(sandboxDir,'eyeTrackingDEMO.zip');
if ~exist (demoPackage,'file')
    url = 'https://ndownloader.figshare.com/files/9355459?private_link=011191afe46841d2c2f5';
    system (['curl -L ' sprintf(url) ' > ' sprintf(demoPackage)])
    currentDir = pwd;
    cd (sandboxDir)
    unzip(demoPackage)
    cd (currentDir)
end


%% Prepare analysis parameters
% We can conceptually break the analysis into two components:
% - Initital processing to the stage of fitting an ellipse to the pupil
%   perimeter
% - Definition of the scene geometry and re-fitting of the pupil
%   perimeter using scene constraints and empirical Bayesian temporal
%   smoothing.


%% Initial processing -- deinterlace to minimally constrained ellipse fitting
% These parameters are used in the function findPupilPerimeter; definition
% of the parameters may be found in the header comments for that routine.
% To select good parameters for your video, use the interactive routine
% estimateSceneParamsGUI.m, which is found in the Utilities directory
pupilFrameMask = [64 109 75 183];
glintFrameMask = [157 148 173 192];
pupilRange = [34 51];
pupilCircleThresh = 0.0179;
pupilGammaCorrection = 0.75;


%% Run the analysis pipeline
% This call to the pipeline will produce the initial ellipse fit and a fit
% video. It includes the stages:
%   deinterlaceVideo
%   findGlint
%   findPupilPerimeter
%   makeControlFile
%   applyControlFile
%   fitPupilPerimeter
%   makeFitVideo
%
runVideoPipeline( pathParams, ...
    'nFrames',nFrames,'verbose', verbose, 'tbSnapshot',tbSnapshot, 'useParallel',true, ...
    'pupilFrameMask', [64 109 75 183], 'glintFrameMask', [157 148 173 192], ...
    'pupilRange', [34 51], 'pupilCircleThresh', 0.0179, 'pupilGammaCorrection', 0.75, ...
    'overwriteControlFile', true, 'catchErrors', false,...
    'lastStageByNumber',6,'makeFitVideoByNumber',6);

% Note that each stage could be called separately, instead of using the
% pipeline command:
%{
    deinterlaceVideo(rawVideoFileName, deinterlacedVideoFileName);
    findGlint(deinterlacedVideoFileName, glintFileName,'glintFrameMask', [157 148 173 192]);
    findPupilPerimeter(deinterlacedVideoFileName, perimeterFileName, ...
        'pupilFrameMask', [64 109 75 183], 'glintFrameMask', [157 148 173 192], ...
        'pupilRange', [34 51], 'pupilCircleThresh', 0.0179, 'pupilGammaCorrection', 0.75);
    makeControlFile(controlFileName,perimeterFileName,glintFileName,...
        'useParallel',true);
    applyControlFile(perimeterFileName,controlFileName,correctedPerimeterFileName);
    fitPupilPerimeter(correctedPerimeterFileName,pupilFileName,...
        'useParallel',true);
    makeFitVideo(deinterlacedVideoFileName,fitVideoFileName,...
        'glintFileName',glintFileName,...
        'controlFileName',controlFileName,...
        'perimeterFileName',correctedPerimeterFileName,...
        'pupilFileName',pupilFileName)
%}


%% Secondary processing -- scene definition and constrained fitting
% Improved fitting can be obtained by defining properties of the eye, the
% camera, and their relationship. Parameters that define the camera and the
% eye described in greater detail in the routine createSceneGeometry which
% is found in the gkaModelEye (https://github.com/gkaguirrelab/gkaModelEye)
% code repository.

% Define camera parameters. These were obtained by an empirical measurement
% (camera resectioning) of the IR camera used to record the demo data. Use
% the matlab routine cameraCalibrator:
%	https://www.mathworks.com/help/vision/ug/single-camera-calibrator-app.html
intrinsicCameraMatrix = [2627.0 0 338.1; 0 2628.1 246.2; 0 0 1];
radialDistortionVector = [-0.3517 3.5353];
sensorResolution = [640 480];
spectralDomain = 'nir';

% Define properties of the eye of the subject. The "maxIrisDiamPixels" is
% the largest horizontal diameter of the iris that is seen in the video.
% This valye is used to estimate the distance of the camera from the eye.
% The spherical ametropia is the refractive error correction (in diopters)
% of the eye of the subject. A negative value is the correction for a
% myopic (near-sighted) person. The axial length would be obtained from an
% ophthalmologic device such as the IOL Master, as would the
% measuredCornealCurvature. If any of these last three values are
% unavailable, leave these parameters undefined.
eyeLaterality = 'right';
maxIrisDiamPixels = 270;
sphericalAmetropia = -1.5;
axialLength = 25.35;
measuredCornealCurvature = [41.36,41.67,25];

% Estimate camera distance from the visible iris diameter in pixels.
% Because biological variation in the size of the visible iris is known, we
% can use the observed iris size in pixels to obtain a guess as to the
% distance of the eye from the camera in mm.
sceneGeometry = createSceneGeometry(...
    'radialDistortionVector',radialDistortionVector, ...
    'intrinsicCameraMatrix',intrinsicCameraMatrix);
[cameraDepthMean, cameraDepthSD] = depthFromIrisDiameter( sceneGeometry, maxIrisDiamPixels );

% Assemble the scene parameter bounds. These are in the order of:
%   [torsion; x; y; z; eyeRotationScalarJoint; eyeRotationScalerDifferential]
% where torsion specifies the torsion of the camera with respect to the eye
% in degrees, [x y z] is the translation of the camera w.r.t. the eye in
% mm, and the eyeRotationScalar variables are multipliers that act upon the
% centers of rotation estimated for the eye. If the eye is markedly
% off-center in the image, then the translation bounds should be increased.
% The upper and lower bounds are fixed here so that the demo completes
% rapidly.
sceneParamsLB =  [3.60; -0.10; -0.70; 109.40; 0.84; 0.96];
sceneParamsLBp = [3.60; -0.10; -0.70; 109.40; 0.84; 0.96];
sceneParamsUBp = [3.60; -0.10; -0.70; 109.40; 0.84; 0.96];
sceneParamsUB =  [3.60; -0.10; -0.70; 109.40; 0.84; 0.96];

% To select good scene parameters for your video, use the interactive
% routine estimateSceneParamsGUI.m, which is found in the Utilities
% directory. You will first have to create a sceneGeometry file.
%{
    sceneGeometryFileName = fullfile(pathParams.dataOutputDirFull,[pathParams.runName '_sceneGeometry.mat']);
    ellipseArrayList = [ 732, 896, 1023, 1167, 1261, 1383, 1542, 1646, 1808 ];
    estimateSceneParamsGUI(sceneGeometryFileName,'ellipseArrayList',ellipseArrayList)
%}


% The estimation of scene geometry is greatly aided by having the subject
% fixate targets at known visual angle positions. The routine
% estimateSceneParams accepts a list of frames of the video during which
% the pupil is well visualized and a list of target positions (X,Y) in
% units of degrees corresponding to each of those video frames.
fixationTargetArray = [ -7, 0, 7, 0, 0, 7, -7, 7, -7 ; 7, -7, 0, 0, 7, 7, 0, -7, -7];
ellipseArrayList = [ 732, 896, 1023, 1167, 1261, 1383, 1542, 1646, 1808 ];

% The routine estimateSceneGeometry performs multiple non-linear searches
% from different x0 start points. Each search is time consuming (e.g., 30
% minutes). They will run in parallel if multiple cores are available. We
% will limit then number of searches to only 2 in this demo, but closer to
% 10 searches yields improved results.
nBADSsearches = 2;


%% Run the analysis pipeline
% Run the video pipeline from the stage of estimation of scene geometry
% through to the end. It includes the stages:
%   estimateSceneParams
%   fitPupilPerimeter (now with scene geometry constraints)
%   smoothPupilRadius
%   makeFitVideo
%
runVideoPipeline( pathParams, ...
    'nFrames',nFrames,'verbose', verbose, 'tbSnapshot',tbSnapshot, 'useParallel',true, ...
    'intrinsicCameraMatrix',intrinsicCameraMatrix, ...
    'radialDistortionVector',radialDistortionVector, ...
    'sensorResolution',sensorResolution,...
    'spectralDomain',spectralDomain, ...
    'eyeLaterality',eyeLaterality,'sphericalAmetropia',sphericalAmetropia,...
    'axialLength',axialLength,'measuredCornealCurvature',measuredCornealCurvature,...
    'sceneParamsLB',sceneParamsLB,'sceneParamsUB',sceneParamsUB,...
    'sceneParamsLBp',sceneParamsLBp,'sceneParamsUBp',sceneParamsUBp,...
    'nBADSsearches',nBADSsearches,...
    'catchErrors', false,...
    'fixationTargetArray',fixationTargetArray,'ellipseArrayList', ellipseArrayList, ...
    'skipStageByNumber',1:6);




%% Plot some fits
pupilFileName = fullfile(pathParams.dataOutputDirFull,[pathParams.runName '_pupil.mat']);
load(pupilFileName,'pupilData');

pupilFileName = fullfile(pathParams.dataOutputDirFull,[pathParams.runName '_timebase.mat']);
load(pupilFileName,'timebase');

% Convert timebase from msecs to minutes
timebase.values = timebase.values./1000./60;


% Identify the time points with good fitting results
rmseThreshold = 3;
highRMSE = pupilData.radiusSmoothed.ellipses.RMSE > rmseThreshold;
fitAtBound = pupilData.radiusSmoothed.eyePoses.fitAtBound;
good = logical(~highRMSE .* ~fitAtBound);

yRangeIncrement = [5 5 0.25];
xLim = [0 ceil(max(timebase.values)*10)/10];

% Set up a figure
figure
eyePoseParamsToPlot = [1 2 4];
yAxisLabels={'azimuth [deg]','elevation [deg]','stop radius [mm]'};

% Loop over the three eyePose parameters to be plotted
for kk=1:length(eyePoseParamsToPlot)
    
    % Get the y limits for this parameter
    lb = floor(min(pupilData.radiusSmoothed.eyePoses.values(good,eyePoseParamsToPlot(kk))) ./ yRangeIncrement(kk)).*yRangeIncrement(kk);
    ub = ceil(max(pupilData.radiusSmoothed.eyePoses.values(good,eyePoseParamsToPlot(kk))) ./ yRangeIncrement(kk)).*yRangeIncrement(kk);
    
    % Define the subplot for this acqusition
    subplot(3,1,kk,'align');
    
    % Plot the time-series. Make the red fit dots transparent
    plot(timebase.values,pupilData.sceneConstrained.eyePoses.values(:,eyePoseParamsToPlot(kk)),'-','Color',[0.85 0.85 0.85],'LineWidth',1);
    hold on
    hLineRed = plot(timebase.values(good),pupilData.radiusSmoothed.eyePoses.values(good,eyePoseParamsToPlot(kk)),'o','MarkerSize',3);
    drawnow
    hMarkerRed = hLineRed.MarkerHandle;
    hMarkerRed.FaceColorData = uint8(255*[1; 0; 0; 0.25]);
    hMarkerRed.FaceColorType = 'truecoloralpha';
    hMarkerRed.EdgeColorData = uint8([0; 0; 0; 0]);
    
    % Add the markers for high RMSE plot points
    lowY = lb + (ub-lb)/20;
    hLineGray = plot(timebase.values(highRMSE),repmat(lowY,size(timebase.values(highRMSE))),'o','MarkerSize',3);
    drawnow
    if ~isempty(hLineGray)
        hMarkerGray = hLineGray.MarkerHandle;
        hMarkerGray.FaceColorData = uint8(255*[0; 0.75; 0; .5]);
        hMarkerGray.FaceColorType = 'truecoloralpha';
        hMarkerGray.EdgeColorData = uint8([0; 0; 0; 0]);
    end
    
    % Add the markers for at bound plot points
    if isfield(pupilData.radiusSmoothed.eyePoses,'fitAtBound')
        hLineBlue = plot(timebase.values(fitAtBound),repmat(lowY,size(timebase.values(fitAtBound))),'o','MarkerSize',3);
        drawnow
        if ~isempty(hLineBlue)
            hMarkerBlue = hLineBlue.MarkerHandle;
            hMarkerBlue.FaceColorData = uint8(255*[0; 0; 0.75; 1]);
            hMarkerBlue.FaceColorType = 'truecoloralpha';
            hMarkerBlue.EdgeColorData = uint8([0; 0; 0; 0]);
        end
    end
    
    % Set the plot limits
    xlim(xLim);
    xticks(fix(xLim(1)):0.25:ceil(xLim(2)*4)/4)
    ylim([lb ub]);
       
    % Add a y-axis label
    ylabel(yAxisLabels{kk});
    
    %   set(gca,'TickDir','out')
    if kk == 1
        title('Demo data transparentTrack');
    end
    if kk ~= length(eyePoseParamsToPlot)
        set(gca,'XColor','none')
    else
        xlabel('time from scan start [minutes]');
    end
    
    box off
    
end % loop over eyePose params


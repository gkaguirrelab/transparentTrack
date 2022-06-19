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
% We can conceptually break the analysis into two components:
% - Initital processing to the stage of fitting an ellipse to the pupil
%   perimeter
% - Estimation of the sceneGeometry
% - re-fitting of the pupil perimeter using scene constraints and empirical
%   Bayesian temporal smoothing.
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

% Download the data if it is not already there. If this automated download
% approach does not work, you can download the file manually by pasting the
% URL into a browser. Place the zip file into the sandboxDir and unzip
% there.
demoPackage = fullfile(sandboxDir,'eyeTrackingDEMO.zip');
if ~exist (demoPackage,'file')
    url = 'https://ndownloader.figshare.com/files/9355459?private_link=011191afe46841d2c2f5';
    system (['curl -L ' sprintf(url) ' > ' sprintf(demoPackage)])
    currentDir = pwd;
    cd (sandboxDir)
    unzip(demoPackage)
    cd (currentDir)
end


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

% Run the analysis pipeline
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
% Note that each stage could be called separately, instead of using the
% pipeline command
runVideoPipeline( pathParams, ...
    'nFrames',nFrames,'verbose', verbose, 'tbSnapshot',tbSnapshot, 'useParallel',true, ...
    'pupilFrameMask', [64 109 75 183], 'glintFrameMask', [157 148 173 192], ...
    'pupilRange', [34 51], 'pupilCircleThresh', 0.0179, 'pupilGammaCorrection', 0.75, ...
    'overwriteControlFile', true, 'catchErrors', false,...
    'lastStageByNumber',6,'makeFitVideoByNumber',6);

    
%% Estimate a sceneGeometry -- TIME CONSUMING
% Improved fitting can be obtained by defining properties of the eye, the
% camera, and their relationship. Parameters that define the camera and the
% eye described in greater detail in the routine createSceneGeometry which
% is found in the gkaModelEye (https://github.com/gkaguirrelab/gkaModelEye)
% code repository. Here, we search over possible parameters of the eye and
% scene model to best fit the set of video images of the eye posed to gaze
% at different targets.

% Define camera parameters. These were obtained by an empirical measurement
% (camera resectioning) of the IR camera used to record the demo data. Use
% the matlab routine cameraCalibrator:
%	https://www.mathworks.com/help/vision/ug/single-camera-calibrator-app.html
intrinsicCameraMatrix = [2627.0 0 338.1; 0 2628.1 246.2; 0 0 1];
radialDistortionVector = [-0.3517 3.5353];
sensorResolution = [640 480];
spectralDomain = 'nir';

% Define properties of the eye of the subject. The spherical ametropia is
% the refractive error correction (in diopters) of the eye of the subject.
% A negative value is the correction for a myopic (near-sighted) person. If
% either of the last two values are unavailable, leave these parameters
% undefined.
eyeLaterality = 'right';
sphericalAmetropia = -1.5;
axialLength = 25.35;

% INTERACTIVE: Estimate camera distance from the visible iris diameter in
% pixels. Because biological variation in the size of the visible iris is
% known, we can use the observed iris size in pixels to obtain a guess as
% to the distance of the eye from the camera in mm. The value is hard-coded
% here, but you can run the commented function "estimateCameraDepth" to
% obtain this value yourself
grayVideoName = fullfile(pathParams.dataOutputDirFull,[pathParams.runName '_gray.avi']);
cameraDepth = 110.22;
%{
    cameraDepth = estimateCameraDepth( grayVideoName );
%}
cameraTranslation = [0; 0; cameraDepth];

% INTERACTIVE: Estimate camera torsion from medial and lateral canthii. The
% angle between the medial and lateral canthus can guide as to the angle of
% rotation of the camera with respect to the head. The value is hard-coded
% here, but you can run the commented function "estimateCameraTorsion" to
% obtain this value yourself
subjectDemographic = 'Caucasian';
cameraTorsionStruct.([eyeLaterality 'Eye' subjectDemographic]) = 4;
%{
    cameraTorsionStruct = estimateCameraTorsion( grayVideoName );
%}
cameraTorsion = cameraTorsionStruct.([eyeLaterality 'Eye' subjectDemographic]);

% Because we have a good estimate of cameraDepth and cameraTorsion, we tell
% the sceneGeometry search to lock these parameters, as well as the primary
% position of the eye. We also limit the search over parameters of the eye,
% including the corneal curvature and the "common camera depth" parameter.
% If we had multiple gazeCal measurements, we could allow greater
% flexibility for these.
%
% We are mostly allowing for a search over camera translation, the rotation
% center of the eye, and a bit of the corneaAxialRadius.
model.scene.bounds = [0 0 0 20 20 0];
model.eye.bounds = [2, 0, 0, 0, 0, 0, 0.25, 0.25, 0];

% The estimation of scene geometry is greatly aided by having the subject
% fixate targets at known visual angle positions. The routine
% estimateSceneParams accepts a list of frames of the video during which
% the pupil is well visualized and a list of target positions (X,Y) in
% units of degrees corresponding to each of those video frames.
gazeTargets = [ -7, 0, 7, 0, 0, 7, -7, 7, -7 ; 7, -7, 0, 0, 7, 7, 0, -7, -7];
frameSet = [ 732, 896, 1023, 1167, 1261, 1383, 1542, 1646, 1808 ];

% The routine estimateSceneGeometry performs multiple non-linear searches
% from different x0 start points. Each search is time consuming (e.g., 30
% minutes). They will run in parallel if multiple cores are available. We
% will limit then number of searches to only 2 in this demo, but closer to
% 10 searches yields improved results.
nBADSsearches = 2;

% Assemble the eye and scene properties
eyeArgs = {'axialLength',axialLength,'sphericalAmetropia',sphericalAmetropia,'eyeLaterality',eyeLaterality,};
sceneArgs = {'intrinsicCameraMatrix',intrinsicCameraMatrix,'radialDistortionVector',radialDistortionVector,...
    'sensorResolution',sensorResolution,'spectralDomain',spectralDomain};

% Run the search
videoStemName = fullfile(pathParams.dataOutputDirFull,pathParams.runName);
sceneObjects = estimateSceneParams(videoStemName, frameSet, gazeTargets,...
    'cameraDepth',cameraDepth,'cameraTorsion',cameraTorsion,...
    'model',model,...
    'eyeArgs',eyeArgs,'sceneArgs',sceneArgs);


%% Apply sceneGeometry to refine pupil boundaries
% Run the second stage of the video pipeline:
%   fitPupilPerimeter (now with scene geometry constraints)
%   smoothPupilRadius (which is a Bayesian temporal smoothing of pupil
%          diameter)
%   makeFitVideo
%
runVideoPipeline( pathParams, ...
    'nFrames',nFrames,'verbose', verbose, 'tbSnapshot',tbSnapshot, 'useParallel',true, ...
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


%  DEMO rendered eye

%% set paths and make directories
% create test sandbox on desktop
sandboxDir = '~/Desktop/eyeTrackingDEMO/renderedEye';
if ~exist(sandboxDir,'dir')
    mkdir(sandboxDir)
end

%% define some names
renderedVideoName = 'pupil_movie';

grayVideoName = fullfile(sandboxDir, [renderedVideoName '.avi']);
glintFileName = fullfile(sandboxDir, [renderedVideoName '_glint.mat']);
perimeterFileName = fullfile(sandboxDir, [renderedVideoName '_perimeter.mat']);
controlFileName = fullfile(sandboxDir, [renderedVideoName '_controlFile.csv']);
correctedPerimeterFileName = fullfile(sandboxDir, [renderedVideoName '_correctedPerimeter.mat']);
pupilFileName = fullfile(sandboxDir, [renderedVideoName '_pupil.mat']);
sceneGeometryFileName = fullfile(sandboxDir, [renderedVideoName '_sceneGeometry.mat']);
sceneDiagnosticPlotFileName = fullfile(sandboxDir, [renderedVideoName '_sceneDiagnosticPlot.pdf']);
finalFitVideoName = fullfile(sandboxDir, [renderedVideoName '_finalFit.avi']);


%% run the pipleline steps

findGlint(grayVideoName, glintFileName,'verbosity','full')

findPupilPerimeter(grayVideoName, perimeterFileName,'verbosity','full');

makeControlFile(controlFileName, perimeterFileName, glintFileName,'verbosity','full');

applyControlFile(perimeterFileName,controlFileName,correctedPerimeterFileName,'verbosity','full');

fitPupilPerimeter(correctedPerimeterFileName, pupilFileName,'nSplits', 0,'verbosity','full');

estimateSceneGeometry(pupilFileName, sceneGeometryFileName, ...
    'sceneDiagnosticPlotFileName', sceneDiagnosticPlotFileName,'verbosity','full');

makeControlFile(controlFileName, perimeterFileName, glintFileName, ...
            'sceneGeometryFileName', sceneGeometryFileName, 'overwriteControlFile', true,'verbosity','full');
            
applyControlFile(perimeterFileName,controlFileName,correctedPerimeterFileName,'verbosity','full');

fitPupilPerimeter(correctedPerimeterFileName, pupilFileName, ...
            'sceneGeometryFileName', sceneGeometryFileName,'verbosity','full');
        
smoothPupilArea(correctedPerimeterFileName, pupilFileName, sceneGeometryFileName,'verbosity','full');

makeFitVideo(grayVideoName, finalFitVideoName, ...
            'glintFileName', glintFileName, 'perimeterFileName', correctedPerimeterFileName,...
            'pupilFileName', pupilFileName, 'sceneGeometryFileName', sceneGeometryFileName, ...
            'controlFileName',controlFileName,'verbosity','full');

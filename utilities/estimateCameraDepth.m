function [cameraDepthMean, cameraDepthSD] = estimateCameraDepth( grayVideoName, varargin )
% Estimate camera depth given sceneGeometry and iris diameter in pixels
%
% Syntax:
%  [cameraDepthMean, cameraDepthSD] = estimateCameraDepth( grayVideoName )
%
% Description:
%   There is limited individual variation in the horizontal visible
%   diameter of the human iris. The maximum observed diameter of the border
%   of the iris in a set of images of the eye will correspond to the
%   diameter of the iris when the eye is posed so a line that connects the
%   center of rotation of the eye and the center of the iris is aligned
%   with the optical axis of the camera. This routine takes a video file
%   and has the user adjust a line to correspond to the maximal diameter of
%   the iris.
%
% Inputs:
% Inputs:
%	grayVideoName         - Full path to the video in which to track the
%
% Optional key/value pairs (flow control)
%  'startFrame'           - First frame from which to start the analysis.
%
% Outputs;
%   cameraDepthMean       - Scalar. The distance, in mm, of the camera
%                           from the corneal apex at pose angle [0, 0, 0],
%                           calculated assuming the average horizontal
%                           visible iris diameter.
%   cameraDepthSD         - Scalar. The 1SD variation in distance, in mm, 
%                           of the camera from the corneal apex at pose
%                           angle [0, 0, 0], given the variation in iris
%                           diameters observed in a population.
%


%% parse input and define variables

p = inputParser; p.KeepUnmatched = true;

% Required
p.addRequired('grayVideoName',@isstr);

% Optional flow control params
p.addParameter('startFrame',1,@isnumeric);

% parse
p.parse(grayVideoName, varargin{:})


% If the grayVideoName is empty, provide a file picker GUI
if isempty(grayVideoName)
    [fileName, path] = uigetfile({'*.mp4;*.mov;*avi'});
    if isempty(fileName)
        returne
    end
    grayVideoName = [path, fileName];
else
    grayVideoName = p.Results.grayVideoName;
end

% Prepare the video object
videoInObj = videoIOWrapper(grayVideoName,'ioAction','read');

% Grab the image dimensions
width = videoInObj.Width;
height = videoInObj.Height;

% read the video desired frame into memory, adjust gamma
thisFrame = read(videoInObj,p.Results.startFrame);

% close the video object
clear videoInObj

% create a figure for display
figureHandle=figure();

% Show the image
imshow(thisFrame,'Border', 'tight', 'InitialMagnification', 200);

roi = images.roi.Line(gca,'Position',[round(width*0.8), round(height*0.25);round(width*0.25),round(height*0.25)]);

% Enter a while loop
notDoneFlag = true;

% Wait until the user is done (presses return or hits 'q')
while notDoneFlag
    
    % Wait for operator input
    waitforbuttonpress
    keyChoiceValue = double(get(gcf,'CurrentCharacter'));
    if ~isempty(keyChoiceValue)
        switch keyChoiceValue
            case {13,27}
                notDoneFlag = false;
        end
    end
end

% The roi position is given as Xpos, Ypos, width, height. [Xpos, Ypos] of
% [0, 0] is the upper left of the displayed image.
x = roi.Position;

% Get the observed iris diameter in pixels
observedIrisDiamPixels = norm(x(1,:)-x(2,:));

% Clean up
close(figureHandle)

% These are the assume sizes of the width of the iris in mm, derived from
% the example code below
trueIrisSizeMean = 5.55;
trueIrisSizeSD = 0.33;

% We now identify the camera distances corresponding the mean, and then to
% the +1SD iris sizes

% Set the x0 position for the search to be the passed scene geometry
sceneGeometry = createSceneGeometry();
x0 = sceneGeometry.cameraPosition.translation(3);

searchSDvals = [0, 1];
for tt = 1:2
    assumedIrisRadius = trueIrisSizeMean + trueIrisSizeSD*searchSDvals(tt);
    cameraTranslationValues(tt) = fminsearch(@objfun, x0);
end
    function fVal = objfun(x)
        candidateSceneGeometry = sceneGeometry;
        candidateSceneGeometry.eye.iris.radius = assumedIrisRadius;
        candidateSceneGeometry.cameraPosition.translation(3) = x;
        [~, ~, imagePoints, ~, ~, ~, pointLabels] = ...
            projectModelEye([0 0 0 1], candidateSceneGeometry, 'fullEyeModelFlag', true, 'nIrisPerimPoints', 20);
        idx = find(strcmp(pointLabels,'irisPerimeter'));
        predictedIrisDiamPixels = max(imagePoints(idx,1))-min(imagePoints(idx,1));
        fVal = (predictedIrisDiamPixels - observedIrisDiamPixels)^2;
    end

cameraDepthMean = cameraTranslationValues(1);
cameraDepthSD = cameraTranslationValues(2)-cameraTranslationValues(1);

end




%% Iris width values
% Define the iris radius. One study measured the horizontal visible
% iris diameter (HVID) in 200 people, and found a mean of 11.8 with
% a range of 10.2 - 13.0.
%    PJ Caroline & MP Andrew. "The Effect of Corneal Diameter on
%    Soft Lens Fitting, Part 1" Contact Lens Spectrum, Issue: April
%    2002
%    https://www.clspectrum.com/issues/2002/april-2002/contact-lens-case-reports
%
% Bernd Bruckner of the company Appenzeller Kontaktlinsen AG
% supplied me with a tech report from his company (HVID & SimK
% study) that measured HVID in 461 people. These data yield a mean
% iris radius of 5.92 mm, 0.28 SD. The values from the histogram
% are represented here, along with a Gaussian fit to the
% distribution
%{
	counts = [0 2 2 0 0 4 5 12 19 23 36 44 52 41 39 37 43 30 28 12 15 10 4 1 0 2 0];
	HVIDRadiusmm = (10.5:0.1:13.1)/2;
	hvidGaussFit = fit(HVIDRadiusmm', counts', 'gauss1');
	hvidRadiusMean = hvidGaussFit.b1;
    hvidRadiusSD =  hvidGaussFit.c1;
    figure
    plot(HVIDRadiusmm, hvidGaussFit(HVIDRadiusmm), '-r')
    hold on
    plot(HVIDRadiusmm, counts, '*k')
    xlabel('HVID radius in mm')
    ylabel('counts')
%}
% The HVID is the refracted iris size. We can use the forward model
% to find the size of the true iris for the mean and mean+1SD observed
% refracted iris sizes
%{
    % Mean and SD radius value from prior block of code
    hvidRadiusMean = 5.9161;
    hvidRadiusSD = 0.2830;
    sceneGeometry = createSceneGeometry();
    % Inactivate ray tracing
    sceneGeometry.refraction = [];
    % Get the area in pixels of a "pupil" that is the same radius
    % as the mean HVID when there is no ray tracing
    hvidP=projectModelEye([0 0 0 hvidRadiusMean],sceneGeometry);
    % Restore ray tracing
    sceneGeometry = createSceneGeometry();
    % Set up the objective function
    myArea = @(p) p(3);
    myObj = @(r) (hvidP(3) - myArea(projectModelEye([0 0 0 r],sceneGeometry)))^2;
    [r,pixelError] = fminsearch(myObj,5.5);
    fprintf('An unrefracted iris radius of %4.2f yields a refracted HVID of %4.2f \n',r,hvidRadiusMean)
    % Now handle the +1SD case
    % Inactivate ray tracing
    sceneGeometry.refraction = [];
    % Get the area in pixels of a "pupil" that is the same radius
    % as the mean HVID when there is no ray tracing
    hvidP=projectModelEye([0 0 0 hvidRadiusMean+hvidRadiusSD],sceneGeometry);
    % Restore ray tracing
    sceneGeometry = createSceneGeometry();
    % Set up the objective function
    myArea = @(p) p(3);
    myObj = @(r) (hvidP(3) - myArea(projectModelEye([0 0 0 r],sceneGeometry)))^2;
    [r,pixelError] = fminsearch(myObj,5.5);
    fprintf('An unrefracted iris radius of %4.2f yields a refracted HVID of %4.2f \n',r,hvidRadiusMean+hvidRadiusSD)
    
%}

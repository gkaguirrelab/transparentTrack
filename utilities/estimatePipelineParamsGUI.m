function [ initialParams ] = estimatePipelineParamsGUI(grayVideoName, approach, varargin)
% Extract parameters to guide pupil tracking by asking  operator input
%
% Syntax:
%  [ initialParams ] = estimatePipelineParamsGUI(grayVideoName, protocol)
%
% Description:
%   Successful pupil tracking requires the input of various parameters,
%   some of which vary across subjects. This routine is intended as a quick
%   means of figuring out what those parameters that vary should be.
%
%   Specifically we are interested in determining: 1) pupilFrameMask
%   (defines the window in which the tracking routine will look for the
%   pupil), 2) glintFrameMask (defines the window in which the tracking
%   routine will look for the glint/glints), 3) pupilRange (sets upper and
%   lower limits for the radius of the pupil the routine searches for in
%   the first frame), 4) pupilCircleThresh (a value which is used to
%   differentiate which intensities correspond to pupil vs. surrounding
%   iris), and 5) maximumVisibleIrisDiameter (defines the diameter of the
%   iris seen in the frame in which the eye is looking straight ahead).
%
%   The process of this routine is as follows. First, the inputted video is
%   opened (using the system's default video program) allowing the user to
%   assess the quality of the video as well as find a representative frame
%   from which we can extract our parameters. The default frame shown first
%   is the first frame, but the operator is prompted to specify if a
%   different frame is required. After the frame has been selected, the
%   operator is instructed to click repeatedly along the boundary of the
%   pupil in the figure window to define the shape of the pupil; a green
%   ellipse fitted to those points is then presented. Next, the operator is
%   asked to click in the figure to define the position of the glints.
%   Unless the routine is changed from its default behavior, the operator
%   is then shown a histogram showing the intensities of pixels found in
%   both the pupil and surrounding iris. On this histogram, the user is
%   asked to select an intensity value that best differetiates iris from
%   pupil. Finally, the user is asked to select a frame in which the
%   recorded eye is looking straight ahead, and the iris therefore is its
%   maximal apparent size. Once this frame has been selected, the operator
%   must click twice on the outer limit of the iris to define its diameter.
%
%   With this process completed, the pupil perimeter found in four frames
%   from the video in question will be displayed, allowing the user to
%   assess the quality of future tracking. Assuming these frames look
%   reasonable, the hope is the default parameters for a given protocol,
%   combined with these newly estimated parameters, should provide
%   reasonable tracking.
%
% Inputs:
%   grayVideoName         - String vector. The full path to a video file
%                           for which parameters are to be derived. If an
%                           empty string is passed (''), the operator can
%                           choose the relevant video file (mp4 or mov) via
%                           a dialog box.
%   approach              - String vector. Defines the default param values
%                           for a particular protocol. Defined values are:
%                           {'TOME','SquintToPulse'}
%
% Optional key/value pairs that control general tracking:
%  'frameNumber'          - A number. Controls which frame of the video is
%                           presented, and upon which the operator will
%                           specify the location of the pupil and glints.
%                           The default value is 1.
%  'openVideo'            - A logical. If true, the relevant video will be
%                           opened by the default video player of the
%                           operator's system.
%  'verbose'              - A logical. Controls the amount of information
%                           displayed onto the console.
%  'ellipseTransparentLB/UB' - A 4-element vector that defines the hard
%                           lower and upper boundaries, respectively, of
%                           the ellipse fitted to the user selected
%                           datapoints around the pupil. The 4 elements, in
%                           order, are: center x position, center y
%                           position, ellipse area, eccentricity, and angle
%                           of tilt.
%  'pupilGammaCorrection' - A number. Gamma correction to be applied to the
%                           video frames (default 1, typical values between
%                           0.5 and 1.8)
%  'frameMaskValue'       - A number. The image value that is assigned to
%                           the region that is masked by frameMask. This
%                           should be a gray that is neither pupil nor
%                           glint.
%  'numberOfGlints'       - A number. Describes how many glints are
%                           expected in this video
%  'maskBox'              - This is the proportion to dilate the pupil
%                           masked region in the vertical and horizontal
%                           directions respectively. A value of 1 will
%                           result in no dilation in that direction. A
%                           value of 2 will result in a masked region that
%                           is twice the size of the pupil radius.
%
% Optional key/value pairs that control estimation of pipeline parameters:
%  'pupilMaskShrinkFactor' - A number. When trying to identify pixels that
%                           are definitively within the pupil, the routine
%                           takes the ellipse fitted to the identified
%                           points along the pupil perimeter and draws a
%                           slightly smaller circle at its center. This
%                           quantity is multiplied by the smaller radius of
%                           the ellipse fit to create the smaller circle.
%  'pupilMaskDilationFactor' - A number. When trying trying to create the
%                           pupilFrameMask, we start from the center of the
%                           user-defined pupil, then expand outwards in
%                           each direction by a length defined by the
%                           radius of the pupil multiplied by the
%                           pupilMaskDilationFactor
%  'inner/outerDilationFactor' - A number. When trying to identify pixels
%                           that are definitively within the pupil, we
%                           dilate the ellipse that defines the pupil by
%                           two separate magnitudes, corresponding to the
%                           inner and outer dilationFactors. Pixels in
%                           between these two dilated circles are
%                           considered iris.
%  'potentialThreshValues' - A vector. Defines the range of values for
%                           pupilCircleThresh that the routine will search
%                           over.
%  'intensityDividerComputeMethod' - A string. Options include 'manual'
%                           (the default option), 'mean', and
%                           'irisMaskMinimum' that refer to various methods
%                           of choosing which pixel intensity value best
%                           differentiates pupil from iris (the variable
%                           intensityDivider in the code below). The
%                           'manual' option has the operator choose the
%                           value from a histogram of pixel intensity
%                           values of both the iris and the pupil. The
%                           'mean' method takes the average pixel intensity
%                           across all pixels in both the iris and pupil
%                           masks as the intensityDivider. The
%                           'irisMaskMinimum' method takes the minimum
%                           pixel intensity within the iris mask as the
%                           intensityDivider.
%  'glintMaskPaddingFactor' - A number. Used to define the limits of the
%                           glintFrameMask. Defines how far from the center
%                           of the identified glints, in units of pixels,
%                           to extend the pupil frame mask in all
%                           directions. Note that the routine multiplies
%                           this value by a scalar to extend the
%                           glintFrameMask in the horizontal direction.
%  'intensityDivider'     - A number. The user has the option to manually
%                           specifiy the intensityDivider value, which
%                           represents the pixel intensity that best
%                           differentiates iris from pupil.
%
% Outputs:
%   initialParams         - A structure with at least 5 subfields, with
%                           subfields including pupilFrameMask,
%                           glintFrameMask, pupilCircleThresh, pupilRange,
%                           and maximumVisibleIrisDiameter. These estimated
%                           parameters an then be added to the default
%                           parameters for a given experiment to begin
%                           processing the relevant video through the
%                           transparentTrack pipeline. Additional subfields
%                           will be added if the operator overrides default
%                           behavior for ellipseTransparentUB/LB,
%                           pupilGammaCorrection, frameMaskValue, or
%                           numberOfGlints

% Examples:
%{
    % ETTBSkip -- This is an interactive example.
    initialParams = estimatePipelineParamsGUI('pathToVideoFile/videoFile.extension')
%}
%{
    % ETTBSkip -- This is an interactive example.
    % Estimate parameters for the TOME dataset.
    initialParams = estimatePipelineParamsGUI('','TOME') 
%}



%% Input parser
p = inputParser; p.KeepUnmatched = true;

% Required
p.addOptional('grayVideoName', [], @(x)(isempty(x) || ischar(x)));
p.addOptional('approach', 'SquintToPulse', @(x)(isempty(x) || ischar(x)));

% Optional flow control params
p.addParameter('frameNumber',3,@isnumeric);
p.addParameter('openVideo',true,@islogical);
p.addParameter('verbose',true,@islogical);

% automatic default parameters in case no default parameters are provided
p.addParameter('ellipseTransparentUB', [], @isnumeric);
p.addParameter('ellipseTransparentLB', [], @isnumeric);
p.addParameter('pupilGammaCorrection', [], @isnumeric);
p.addParameter('frameMaskValue', [], @isnumeric);
p.addParameter('numberOfGlints', [], @isnumeric);
p.addParameter('maskBox', [], @isnumeric);
p.addParameter('smallObjThresh', [], @isnumeric);


% parameters that adjust this initial parameter guessing
p.addParameter('pupilMaskShrinkFactor', 0.9, @isnumeric);
p.addParameter('pupilMaskDilationFactor', 4, @isnumeric);
p.addParameter('pupilRangeDilator', 1.1, @isnumeric);
p.addParameter('pupilRangeContractor', 0.9, @isnumeric);
p.addParameter('innerDilationFactor', 1.1, @isnumeric);
p.addParameter('outerDilationFactor', 1.3, @isnumeric);
p.addParameter('potentialThreshValues', [0.001:0.001:0.2], @isnumeric);
p.addParameter('intensityDividerComputeMethod', 'manual', @isstr);
p.addParameter('glintMaskPaddingFactor', 75, @isnumeric);
p.addParameter('intensityDivider', [], @isnumeric);


% parse
p.parse(grayVideoName,approach, varargin{:})

if strcmp(p.Results.approach, 'TOME')
    ellipseTransparentUB = [1280, 720, 90000, 0.6, pi];
    ellipseTransparentLB = [0, 0, 1000, 0, 0];
    pupilGammaCorrection = 0.75;
    frameMaskValue = 220;
    numberOfGlints = 1;
    maskBox = [2 2];
    smallObjThresh = 400;
elseif strcmp(p.Results.approach, 'SquintToPulse')
    ellipseTransparentUB = [1280, 720, 90000, 0.6, pi];
    ellipseTransparentLB = [0, 0, 1000, 0, 0];
    pupilGammaCorrection = 0.75;
    frameMaskValue = 220;
    numberOfGlints = 2;
    maskBox = [1 1];
    smallObjThresh = 2500;
elseif strcmp(p.Results.approach, 'mtrpGlare')
    ellipseTransparentUB = [640, 480, 90000, 0.8, pi];
    ellipseTransparentLB = [0, 0, 1000, 0, 0];
    pupilGammaCorrection = 0.75;
    frameMaskValue = 220;
    numberOfGlints = 1;
    maskBox = [1 1];
    smallObjThresh = 2500;
elseif strcmp(p.Results.approach, 'SminusMel')
    ellipseTransparentUB = [640, 480, 1e7, 0.99, pi];
    ellipseTransparentLB = [0, 0, 1000, 0, 0];
    pupilGammaCorrection = 0.75;
    frameMaskValue = 220;
    numberOfGlints = 1;
    maskBox = [1 1];
    smallObjThresh = 2500;
    
end

% allow ability to override defaultParams if necessary by passing key-value
% pairs
if ~isempty(p.Results.ellipseTransparentUB)
    ellipseTransparentUB = p.Results.ellipseTransparentUB;
    initialParams.ellipseTransparentUB = p.Results.ellipseTransparentUB;
end
if ~isempty(p.Results.ellipseTransparentLB)
    ellipseTransparentLB = p.Results.ellipseTransparentLB;
    initialParams.ellipseTransparentLB = p.Results.ellipseTransparentLB;
    
end
if ~isempty(p.Results.frameMaskValue)
    frameMaskValue = p.Results.frameMaskValue;
    initialParams.frameMaskValue = p.Results.frameMaskValue;
    
end
if ~isempty(p.Results.pupilGammaCorrection)
    pupilGammaCorrection = p.Results.pupilGammaCorrection;
    initialParams.pupilGammaCorrection = p.Results.pupilGammaCorrection;
    
end
if ~isempty(p.Results.numberOfGlints)
    numberOfGlints = p.Results.numberOfGlints;
    initialParams.numberOfGlints = p.Results.numberOfGlints;
    
end
if ~isempty(p.Results.maskBox)
    maskBox = p.Results.maskBox;
    initialParams.maskBox = p.Results.maskBox;
    
end

if isempty(grayVideoName)
    [fileName, path] = uigetfile({'*.mp4;*.mov;*.avi'});
    grayVideoName = [path, fileName];
end

%% Load up the frame of interest
close all
videoInObj = VideoReader(grayVideoName);
nFrames = floor(videoInObj.Duration*videoInObj.FrameRate);
videoSizeX = videoInObj.Width;
videoSizeY = videoInObj.Height;

% load up the single desired video frame
videoInObj.CurrentTime = (p.Results.frameNumber - 1)/(videoInObj.FrameRate);
thisFrame = readFrame(videoInObj);
thisFrame = rgb2gray(thisFrame);
thisFrame = squeeze(thisFrame);

% open video file, if asked
if p.Results.openVideo
    [recordedErrorFlag, consoleOutput] = system(['open ''' grayVideoName '''']);
end

% present the video frame

figure;
imshow(thisFrame, 'Border', 'tight')
hold on

% ask the user if we the presented frame is OK, or if we should look for a
% different one
if p.Results.frameNumber ==3
    frameCheckChoice = GetWithDefault('>> Is this a good frame? Enter ''y'' to proceed, or ''n'' to choose new frame. [y/n]', 'y');
    if strcmp(frameCheckChoice, 'n')
        close all
        frameRequest = GetWithDefault('>> Enter desired frame:', 1);
        videoInObj.CurrentTime = (frameRequest - 1)/(videoInObj.FrameRate);
        thisFrame = readFrame(videoInObj);
        thisFrame = rgb2gray(thisFrame);
        thisFrame = squeeze(thisFrame);
        figure;
        imshow(thisFrame, 'Border', 'tight');
        hold on
    end
end






%% Guess initial pupil position
% begin user input
fprintf('Define pupil boundary in figure.\n')
string = 'Define the pupil boundary by clicking along the pupil perimeter. Press ENTER when finished.';
hText = text(1,10,string, 'FontSize', 16, 'BackgroundColor', 'white');
[x,y] = ginput;

delete(hText)
% fit an ellipse to the inputted points
[ellipseFitParams] = constrainedEllipseFit(x,y, ellipseTransparentLB, ellipseTransparentUB, []);

% convert the ellipse params from transparent params to explicit params
explicitEllipseFitParams = ellipse_transparent2ex(ellipseFitParams);

% convert the explicit ellipse params to implicit
pFitImplicit = ellipse_ex2im(explicitEllipseFitParams);
% write the implicit function on teh basis of these params
fh=@(x,y) pFitImplicit(1).*x.^2 +pFitImplicit(2).*x.*y +pFitImplicit(3).*y.^2 +pFitImplicit(4).*x +pFitImplicit(5).*y +pFitImplicit(6);
% superimpose the ellipse using fimplicit
hold on
fHandle = fimplicit(fh,[1, videoSizeX, 1, videoSizeY],'Color', 'green','LineWidth',1);
set(gca,'position',[0 0 1 1],'units','normalized')
axis off;

% ask if we need to redo ellipse placement
pupilCheckChoice = GetWithDefault('>> Satisfied with ellipse fit? Enter ''y'' to proceed, or ''n'' to redraw. [y/n]', 'y');
if ~strcmp(pupilCheckChoice, 'y')
    delete(fHandle);
    
    pupilEllipseDoneFlag = false;
    while ~pupilEllipseDoneFlag
        string = 'Define the pupil boundary by clicking along the pupil perimeter. Press ENTER when finished.';
        hText = text(1,10,string, 'FontSize', 16, 'BackgroundColor', 'white');
        [x,y] = ginput;
        
        delete(hText)
        % fit an ellipse to the inputted points
        [ellipseFitParams] = constrainedEllipseFit(x,y, ellipseTransparentLB, ellipseTransparentUB, []);
        
        % convert the ellipse params from transparent params to explicit params
        explicitEllipseFitParams = ellipse_transparent2ex(ellipseFitParams);
        
        % convert the explicit ellipse params to implicit
        pFitImplicit = ellipse_ex2im(explicitEllipseFitParams);
        % write the implicit function on teh basis of these params
        fh=@(x,y) pFitImplicit(1).*x.^2 +pFitImplicit(2).*x.*y +pFitImplicit(3).*y.^2 +pFitImplicit(4).*x +pFitImplicit(5).*y +pFitImplicit(6);
        % superimpose the ellipse using fimplicit
        hold on
        fHandle = fimplicit(fh,[1, videoSizeX, 1, videoSizeY],'Color', 'green','LineWidth',1);
        set(gca,'position',[0 0 1 1],'units','normalized')
        axis off;
        
        pupilCheckChoice = GetWithDefault('>> Satisfied with ellipse fit? Enter ''y'' to proceed, or ''n'' to redraw. [y/n]', 'y');
        
        switch pupilCheckChoice
            case 'y'
                pupilEllipseDoneFlag = true;
                
            case 'n'
                pupilEllipseDoneFlag = false;
                delete(fHandle);
            otherwise
        end
    end
end

%% Figure out the pupil frame

% figure out smaller ellipse axis
if explicitEllipseFitParams(3) > explicitEllipseFitParams(4)
    circleRadius = explicitEllipseFitParams(4);
else
    circleRadius = explicitEllipseFitParams(3);
end

% we'll take the center of the ellipse, and then extend out along from the
% expanded pupil radius to define the pupilFrameMask
if (p.Results.pupilMaskDilationFactor * circleRadius) < 200
    initialParams.pupilFrameMask = [round(explicitEllipseFitParams(2)-200) ...
        round(videoSizeX - (explicitEllipseFitParams(1)+200)) ...
        round(videoSizeY - (explicitEllipseFitParams(2)+200)) ...
        round(explicitEllipseFitParams(1)-200)];
else
    initialParams.pupilFrameMask = [round(explicitEllipseFitParams(2)-explicitEllipseFitParams(4)*p.Results.pupilMaskDilationFactor) ...
        round(videoSizeX - (explicitEllipseFitParams(1)+explicitEllipseFitParams(3)*p.Results.pupilMaskDilationFactor)) ...
        round(videoSizeY - (explicitEllipseFitParams(2)+explicitEllipseFitParams(4)*p.Results.pupilMaskDilationFactor)) ...
        round(explicitEllipseFitParams(1)-explicitEllipseFitParams(3)*p.Results.pupilMaskDilationFactor)];
end

% make sure that pupilFrameMask is actually within the window
for ii = 1:4
    if initialParams.pupilFrameMask(ii) < 1
        initialParams.pupilFrameMask(ii) = 1;
    end
end


%% Figure out initial pupil range
% make the lower bound a little smaller, and upper bound bigger, than the
% user-inputted pupil ellipse
initialParams.pupilRange = [round(circleRadius*p.Results.pupilRangeContractor) round(circleRadius*p.Results.pupilRangeDilator)];

%% Figure out the glinFrameMask
% Ask the user to define the position of the two glints
fprintf('Define glint position in figure.\n')
string = sprintf('Define the glint position by clicking each of the %d glints.', numberOfGlints);
hText = text(1,10,string, 'FontSize', 16, 'BackgroundColor', 'white');
[x,y] = ginput(2);
gHandle = plot(x ,y, '+', 'Color', 'red');

% ask if we need to redo glint placement
glintCheckChoice = GetWithDefault('>> Satisfied with glint locations? Enter ''y'' to proceed, or ''n'' to redraw. [y/n]', 'y');
if ~strcmp(glintCheckChoice, 'y')
    delete(gHandle);
    glintDoneFlag = false;
    while ~glintDoneFlag
        string = sprintf('Define the glint position by clicking each of the %d glints.', numberOfGlints);
        hText = text(1,10,string, 'FontSize', 16, 'BackgroundColor', 'white');
        [x,y] = ginput(2);
        gHandle = plot(x ,y, '+', 'Color', 'red');
        glintCheckChoice = GetWithDefault('>> Satisfied with glint locations? Enter ''y'' to proceed, or ''n'' to redraw. [y/n]', 'y');
        switch glintCheckChoice
            case 'y'
                glintDoneFlag = true;
            case 'n'
                glintDoneFlag = false;
                delete(gHandle);
        end
    end
end

% assume the middle x coordinate of the glintFrameMask should be the
% average of the two x positions of each glint
glintXPosition = mean(x);

% define upper and lower glint positions
glintYPositionLower = max(y);
glintYPositionUpper = min(y);

% expand beyond these positions to define to the glintFrameMask
initialParams.glintFrameMask = [round(glintYPositionUpper - p.Results.glintMaskPaddingFactor) ...
    round(videoSizeX - (glintXPosition + p.Results.glintMaskPaddingFactor*2)) ...
    round(videoSizeY - (glintYPositionLower + p.Results.glintMaskPaddingFactor)) ...
    round(glintXPosition - p.Results.glintMaskPaddingFactor*2)];
% make sure that glintFrameMask is actually within the window
for ii = 1:4
    if initialParams.glintFrameMask(ii) < 1
        initialParams.glintFrameMask(ii) = 1;
    end
end


close all




%% Figure out the pupilCircleThresh
% define a pupil mask. we're interested in figuring out what values make up
% the user-defined pupil

% instantiate mask
pupilMask = zeros(size(thisFrame));

% add circle to mask, corresponding to the position of the user-defined
% circle. we'll make the circle a slightly smaller radius so we're not
% accidentally getting any iris
pupilMask = insertShape(pupilMask,'FilledCircle',[explicitEllipseFitParams(1) explicitEllipseFitParams(2) circleRadius*p.Results.pupilMaskShrinkFactor],'Color','white');
%binarize the mask
pupilMask = im2bw(pupilMask);
% conver thisFrame so we can do some math with it
%thisFrameRGB = rgb2gray(thisFrame);
thisFrameRGB = thisFrame;

% apply the pupilMask to the image, which gives us just the maskedPupil
maskedPupil = immultiply(thisFrameRGB,pupilMask);
% set all values of 0 (ie not pupil) to NaN, for ease of understanding the
% statistics
maskedPupilNaN=double(maskedPupil);
maskedPupilNaN(maskedPupil == 0) = NaN;

% make irisMask
irisMask = zeros(size(thisFrame));
% we'll define two circles, the inner corresponding to the inner iris-pupil
% boundary, the outer corresponding to some circle in the middle of the
% iris. we find these circles by dilating the original pupilCircle, first
% by a little bit (to make the innerIrisMask), then by a little bit more
% (to create the outerIrisMask)
innerIrisMask = insertShape(irisMask,'FilledCircle',[explicitEllipseFitParams(1) explicitEllipseFitParams(2) circleRadius*p.Results.innerDilationFactor],'Color','white');
outerIrisMask = insertShape(irisMask,'FilledCircle',[explicitEllipseFitParams(1) explicitEllipseFitParams(2) circleRadius*p.Results.outerDilationFactor],'Color','white');
% binarize these masks
innerIrisMask = im2bw(innerIrisMask);
outerIrisMask = im2bw(outerIrisMask);
% take the difference, which should just be the iris
differentialIrisMask = outerIrisMask - innerIrisMask;
differentialIrisMask = im2bw(differentialIrisMask);

% apply this mask to the image
maskedIris = immultiply(thisFrameRGB, differentialIrisMask);
maskedIrisNaN = double(maskedIris);
% set all 0 values to NaN, again for ease of statistical analysis
maskedIrisNaN(maskedIris == 0) = NaN;

% the findPupilPerimeter code. first we have to manipulate the image a
% little bit, to mimic some of the settings applied to the relevant region
% as part of the findPupilPerimeter routine. we want to figure out the
% proper value that will be used there, so we'll set up the environment so
% that it'll be the same
% gamma correct the image
thisFrameMasked = imadjust(thisFrame,[],[],pupilGammaCorrection);

% conver to grayscale
%thisFrameMasked = rgb2gray(thisFrameMasked);

% make the image, according to pupilFrameMask
thisFrameMasked((1:initialParams.pupilFrameMask(1)),:) = frameMaskValue; %top
thisFrameMasked(:, (end - initialParams.pupilFrameMask(2):end)) = frameMaskValue; %left
thisFrameMasked((end - initialParams.pupilFrameMask(3):end),:) = frameMaskValue; %bottom
thisFrameMasked(:, (1:initialParams.pupilFrameMask(4))) = frameMaskValue; %right

% smooth the image
I = thisFrameMasked;
filtSize = round([0.01*min(size(I)) 0.01*min(size(I)) 0.01*min(size(I))]);
padP = padarray(I,[size(I,1)/2 size(I,2)/2], 128);
h = fspecial('gaussian',[filtSize(1) filtSize(2)],filtSize(3));
pI = imfilter(padP,h);
pI = pI(size(I,1)/2+1:size(I,1)/2+size(I,1),size(I,2)/2+1:size(I,2)/2+size(I,2));

% figure out the intensity value that differentiates the pupil from the
% iris
if isempty(p.Results.intensityDivider)
    if strcmp(p.Results.intensityDividerComputeMethod, 'mean')
        % one way to do so is to take the average of the values corresponding
        % to the pupil and iris, such that our dividing point is somewhere in
        % between
        intensityDivider = nanmean([maskedPupilNaN(:); maskedIrisNaN(:)]);
    elseif strcmp(p.Results.intensityDividerComputeMethod, 'irisMaskMinimum')
        % alternatively, we can take the dividing point as the minimum value
        % observed in the iris
        intensityDivider = min(maskedIrisNaN(:));
    elseif strcmp(p.Results.intensityDividerComputeMethod, 'manual')
        fprintf('Manually identify intensityDivider in figure window.\n')
        
        plotFig = figure;
        set(gcf,'un','n','pos',[.05,.05,.7,.6])
        axesP = axes('Parent', plotFig);
        hold on
        h1 = histogram(maskedIrisNaN);
        h2 = histogram(maskedPupilNaN);
        
        textYlocation = 1.1*max([h1.Values, h2.Values]);
        legend('Iris', 'Pupil')
        xlabel('Pixel Intensity')
        ylabel('Count')
        string = sprintf('Click once to choose the pixel intensity that best differentiates pupil from iris.');
        hText = text(1,textYlocation, string, 'FontSize', 16, 'BackgroundColor', 'white');
        [x,y] = ginput(1);
        intensityDivider = x;
        close(plotFig);
    else
        warning('Intensity divider compute method not found. Please use either ''mean'', ''manual'', or ''irisMaskMinimum''')
    end
else
    intensityDivider = p.Results.intensityDivider;
end

% i don't know of an analytical way to solve this quantile function, so
% i'll brute-force check what pupilCircleThresh value provides a cutoff
% closest to our intended intensityDivider
counter = 1;
for xx = p.Results.potentialThreshValues
    yy(counter) = quantile(double(pI(:)),xx);
    counter = counter+1;
end
[minValue, potentialIndices] = min(abs(yy-intensityDivider));

% sometimes the attempt to find the potentialThreshValue fails (often
% because the intensityDivider value is not reasonable). here we give the
% user a chance to select the intensityDivider manually, then re-compute
% the potentailThreshValue one last time.
if (min(potentialIndices)-1) < 1 || (min(potentialIndices)-1) > length(p.Results.potentialThreshValues)
    fprintf('Attempt to determine pupilCircleThresh failed. Try again manually\n')
    plotFig = figure;
    set(gcf,'un','n','pos',[.05,.05,.7,.6])
    axesP = axes('Parent', plotFig);
    hold on
    h1 = histogram(maskedIrisNaN);
    h2 = histogram(maskedPupilNaN);
    
    textYlocation = 1.1*max([h1.Values, h2.Values]);
    legend('Iris', 'Pupil')
    xlabel('Pixel Intensity')
    ylabel('Count')
    string = sprintf('Click once to choose the pixel intensity that best differentiates pupil from iris.');
    hText = text(1,textYlocation, string, 'FontSize', 16, 'BackgroundColor', 'white');
    [x,y] = ginput(1);
    intensityDivider = x;
    close(plotFig);
    [minValue, potentialIndices] = min(abs(yy-intensityDivider));
    
end


pupilCircleThresh = p.Results.potentialThreshValues(max([1 min(potentialIndices)-1]));
initialParams.pupilCircleThresh = pupilCircleThresh;

%% Figure out the maximum visible iris diameter
irisDiameterFrame = GetWithDefault('>> Enter the frame number in which the visible iris diameter is largest.', [1]);
videoInObj.CurrentTime = (irisDiameterFrame - 1)/(videoInObj.FrameRate);
thisFrameIris = readFrame(videoInObj);
thisFrameIris = rgb2gray(thisFrameIris);
thisFrameIris = squeeze(thisFrameIris);
figure;
imshow(thisFrameIris, 'Border', 'tight');
hold on
fprintf('Define the maximum visible iris diameter in the figure.\n')
string = sprintf('Define the iris diameter by clicking twice on the outer boundary of the iris.');
hText = text(1,10,string, 'FontSize', 16, 'BackgroundColor', 'white');
[x,y] = ginput(2);
iHandle = plot(x ,y, '+', 'Color', 'red');

% ask if we need to redo glint placement
irisCheckChoice = GetWithDefault('>> Satisfied with iris diameter location? Enter ''y'' to proceed, or ''n'' to redraw. [y/n]', 'y');
if ~strcmp(irisCheckChoice, 'y')
    delete(gHandle);
    irisDoneFlag = false;
    while ~irisDoneFlag
        string = sprintf('Define the iris diameter by clicking twice on the outer boundary of the iris.');
        hText = text(1,10,string, 'FontSize', 16, 'BackgroundColor', 'white');
        [x,y] = ginput(2);
        iHandle = plot(x ,y, '+', 'Color', 'red');
        irisCheckChoice = GetWithDefault('>> Satisfied with iris diameter location? Enter ''y'' to proceed, or ''n'' to redraw. [y/n]', 'y');
        switch irisCheckChoice
            case 'y'
                irisDoneFlag = true;
            case 'n'
                irisDoneFlag = false;
                delete(iHandle);
        end
    end
end
close all

initialParams.maximumVisibleIrisDiameter = abs(x(1) - x(2));

% Could calculate and store the actual default sceneParam values here:
%{
    % Estimate camera distance from iris diameter in pixels
    % Because biological variation in the size of the visible iris is known,
    % we can use the observed maximum diameter of the iris in pixels to obtain
    % a guess as to the distance of the eye from the camera.
    sceneGeometry = createSceneGeometry(...
        'radialDistortionVector',radialDistortionVector, ...
        'intrinsicCameraMatrix',intrinsicCameraMatrix);
    [cameraDepthMean, cameraDepthSD] = depthFromIrisDiameter( sceneGeometry, maxIrisDiamPixels );

    % Assemble the scene parameter bounds. These are in the order of:
    %   torsion, x, y, z, eyeRotationScalarJoint, eyeRotationScalerDifferential
    % where torsion specifies the torsion of the camera with respect to the eye
    % in degrees, [x y z] is the translation of the camera w.r.t. the eye in
    % mm, and the eyeRotationScalar variables are multipliers that act upon the
    % centers of rotation estimated for the eye.
    sceneParamsLB = [-5; -5; -5; cameraDepthMean-cameraDepthSD*2; 0.75; 0.9];
    sceneParamsLBp = [-3; -2; -2; cameraDepthMean-cameraDepthSD*1; 0.85; 0.95];
    sceneParamsUBp = [3; 2; 2; cameraDepthMean+cameraDepthSD*1; 1.15; 1.05];
    sceneParamsUB = [5; 5; 5; cameraDepthMean+cameraDepthSD*2; 1.25; 1.1];
%}




%% diagnostics
% let's see how well we can find the pupil perimeter with these initial
% parameters
nFramesToCheck = 4;
for ii = 1:nFramesToCheck
    framesToCheck(ii) = round((nFrames/(nFramesToCheck-1))*(ii-1))+1;
end
framesToCheck(1) = p.Results.frameNumber;
framesToCheck(end) = nFrames;


counter = 1;
for ii = framesToCheck
    perimeter = [];
%     plotFig = figure;
%     hold on
%     
%     %subplot(2, round(nFrames/2), counter)
%     counter = counter + 1;
%     string = [];
%     string = (['Frame ', num2str(ii)]);
    
    videoInObj.CurrentTime = (ii - 1)/(videoInObj.FrameRate);
    thisFrameDiagnostics = readFrame(videoInObj);
    thisFrameDiagnostics = rgb2gray(thisFrameDiagnostics);
    thisFrameDiagnostics = squeeze(thisFrameDiagnostics);
    
    
    perimeter = findPupilPerimeter(grayVideoName, 'temp', ...
        'startFrame', ii, ...
        'nFrames', 1, ...
        'ellipseTransparentUB', ellipseTransparentUB, ...
        'ellipseTransparentLB', ellipseTransparentLB, ...
        'pupilGammaCorrection', pupilGammaCorrection, ...
        'frameMaskValue', frameMaskValue, ...
        'pupilFrameMask', initialParams.pupilFrameMask, ...
        'pupilRange', initialParams.pupilRange, ...
        'pupilCircleThresh', initialParams.pupilCircleThresh, ...
        'maskBox', maskBox, ...
        'smallObjThresh', smallObjThresh, 'displayMode', true);
%     displayFrame=thisFrameDiagnostics;
%     if ~isempty(perimeter.data{1}.Xp)
%         displayFrame(sub2ind(size(thisFrameDiagnostics),perimeter.data{1}.Yp,perimeter.data{1}.Xp))=255;
%     end
%     if isempty(perimeter.data{1}.Xp)
%         string = 'No pupil found';
%         text(350, 200, string);
%     end
%     imshow(displayFrame, 'Border', 'tight')
%     dText = text(1,10,string, 'FontSize', 16, 'BackgroundColor', 'white');
    delete('temp.mat')
    
    
line([0, size(thisFrameDiagnostics, 2)], [initialParams.glintFrameMask(1), initialParams.glintFrameMask(1)], 'Color', 'r')
line([0, size(thisFrameDiagnostics, 2)], [(size(thisFrameDiagnostics, 1) - initialParams.glintFrameMask(3)), (size(thisFrameDiagnostics, 1) - initialParams.glintFrameMask(3))], 'Color', 'r')
line([initialParams.glintFrameMask(4), initialParams.glintFrameMask(4)], [0, size(thisFrameDiagnostics, 1)], 'Color', 'r');
line([(size(thisFrameDiagnostics, 2) - initialParams.glintFrameMask(2)), (size(thisFrameDiagnostics, 2) - initialParams.glintFrameMask(2))], [0, size(thisFrameDiagnostics, 1)], 'Color', 'r');

end

%% allow the user to adjust certain parameters, then test finding the pupil perimeter again
adjustParamsChoice = GetWithDefault('>> Satisfied with these parameters? Enter ''y'' to proceed and exit, or ''n'' to manually adjust the parameters. [y/n]', 'y');
if ~strcmp(adjustParamsChoice, 'y')
    close all
    adjustParamsFlag = false;
    while ~adjustParamsFlag
        
        fprintf('Select the parameter you would like to adjust:\n')
        fprintf('\t1. ellipseTransparentUB: %g %g %g %g %g \n', ellipseTransparentUB(:));
        fprintf('\t2. ellipseTransparentLB: %g %g %g %g %g \n', ellipseTransparentLB(:));
        fprintf('\t3. pupilGammaCorrection: %g\n', pupilGammaCorrection);
        fprintf('\t4. frameMaskValue: %g\n', frameMaskValue);
        fprintf('\t5. pupilFrameMask: %g %g %g %g\n', initialParams.pupilFrameMask(:));
        fprintf('\t6. pupilCircleThresh: %g\n', initialParams.pupilCircleThresh);
        fprintf('\t7. maskBox: %g %g\n', maskBox(:));
        fprintf('\t8. pupilRange: %g %g\n', initialParams.pupilRange(:));
        fprintf('\t9. smallObjThresh: %g\n', smallObjThresh);

        
        choice = input('\nYour choice: ', 's');
        
        switch choice
            case '1'
                ellipseTransparentUB = input('Enter new ellipseTransparentUB:     ');
                initialParams.ellipseTransparentUB = ellipseTransparentUB;
            case '2'
                ellipseTransparentLB = input('Enter new ellipseTransparentLB:     ');
                initialParams.ellipseTransparentLB = ellipseTransparentLB;
            case '3'
                pupilGammaCorrection = input('Enter new pupilGammaCorrection:     ');
                initialParams.pupilGammaCorrection = pupilGammaCorrection;
            case '4'
                frameMaskValue = input('Enter new frameMaskValue:     ');
                initialParams.frameMaskValue = frameMaskValue;
            case '5'
                initialParams.pupilFrameMask = input('Enter new pupilFrameMask:     ');
            case '6'
                initialParams.pupilCircleThresh = input('Enter new pupilCircleThresh:     ');
            case '7'
                maskBox = input('Enter new maskBox:     ');
                initialParams.maskBox = maskBox;                
            case '8'
                initialParams.pupilRange = input('Enter new pupilRange:     ');
            case '9'
                smallObjThresh = input('Enter new smallObjThresh:       ');
        end
        
        fprintf('New parameters:\n')
        fprintf('\tellipseTransparentUB: %g %g %g %g %g \n', ellipseTransparentUB(:));
        fprintf('\tellipseTransparentLB: %g %g %g %g %g \n', ellipseTransparentLB(:));
        fprintf('\tpupilGammaCorrection: %g\n', pupilGammaCorrection);
        fprintf('\tframeMaskValue: %g\n', frameMaskValue);
        fprintf('\tpupilFrameMask: %g %g %g %g\n', initialParams.pupilFrameMask(:));
        fprintf('\tpupilCircleThresh: %g\n', initialParams.pupilCircleThresh);
        fprintf('\tmaskBox: %g %g\n', maskBox(:));
        fprintf('\tpupilRange: %g %g\n', initialParams.pupilRange(:));

        
        nFramesToCheck = 4;
        for ii = 1:nFramesToCheck
            framesToCheck(ii) = round((nFrames/(nFramesToCheck-1))*(ii-1))+1;
        end
        framesToCheck(1) = p.Results.frameNumber;
        framesToCheck(end) = nFrames;
        
        
        counter = 1;
        for ii = framesToCheck
            perimeter = [];
%             plotFig = figure;
%             hold on
%             
%             %subplot(2, round(nFrames/2), counter)
%             counter = counter + 1;
%             string = [];
%             string = (['Frame ', num2str(ii)]);
            
            videoInObj.CurrentTime = (ii - 1)/(videoInObj.FrameRate);
            thisFrameDiagnostics = readFrame(videoInObj);
            thisFrameDiagnostics = rgb2gray(thisFrameDiagnostics);
            thisFrameDiagnostics = squeeze(thisFrameDiagnostics);
            
            
            perimeter = findPupilPerimeter(grayVideoName, 'temp', ...
                'startFrame', ii, ...
                'nFrames', 1, ...
                'ellipseTransparentUB', ellipseTransparentUB, ...
                'ellipseTransparentLB', ellipseTransparentLB, ...
                'pupilGammaCorrection', pupilGammaCorrection, ...
                'frameMaskValue', frameMaskValue, ...
                'pupilFrameMask', initialParams.pupilFrameMask, ...
                'pupilRange', initialParams.pupilRange, ...
                'pupilCircleThresh', initialParams.pupilCircleThresh, ...
                'maskBox', maskBox, ...
                'smallObjThresh', smallObjThresh,  'displayMode', true);
%             displayFrame=thisFrameDiagnostics;
%             if ~isempty(perimeter.data{1}.Xp)
%                 displayFrame(sub2ind(size(thisFrameDiagnostics),perimeter.data{1}.Yp,perimeter.data{1}.Xp))=255;
%             end
%             if isempty(perimeter.data{1}.Xp)
%                 string = 'No pupil found';
%                 text(350, 200, string);
%             end
%             imshow(displayFrame, 'Border', 'tight')
%             dText = text(1,10,string, 'FontSize', 16, 'BackgroundColor', 'white');
            delete('temp.mat')
        end
        adjustParamsChoice = GetWithDefault('>> Satisfied with these parameters? Enter ''y'' to proceed and exit, or ''n'' to manually adjust the parameters. [y/n]', 'y');
        switch adjustParamsChoice
            case 'y'
                adjustParamsFlag = true;
            case 'n'
                adjustParamsFlag = false;
                close all
        end
        
        
    end
end

%% dump out the params we estimated
if p.Results.verbose
    initialParams
end

end % MAIN


%%% LOCAL FUNCTIONS

function inputVal = GetWithDefault(prompt,defaultVal)
% inputVal = GetWithDefault(prompt,defaultVal)
%
% Prompt for a number or string, with a default returned if user
% hits return.
%
% 4/3/10  dhb  Wrote it.

if (ischar(defaultVal))
    inputVal = input(sprintf([prompt ' [%s]: '],defaultVal),'s');
else
    inputVal = input(sprintf([prompt ' [%g]: '],defaultVal));
end
if (isempty(inputVal))
    inputVal = defaultVal;
end

end % GetWithDefault

% TEST_ellipseThetaConsistency

% This script generates ellipses on the image plane that are
% pseudoPerspective projections of a model pupil. These ellipses are then
% fit with different routines and the resulting output is checked to ensure
% that the convention for specifying the theta value (angle of rotation) of
% the ellipses is the same across routines. Note that this test is
% undefined when the eccentricity of the ellipse is very low, as a circle
% does not have a defined angle of tilt.

%% Clean up
clearvars
close all

% Reset the random number generator so that we have reproducible results
rng default

thetaDifferenceTolerance = 0.05;     % How much variance in theta is allowed?
eccentricityThresh = 0.1;           % Only test ellipses more eccentric than this

% 1. define synthetic data length
nFrames = 200;

% 2. define scene dimension in pixels
videoSizeX = 640;
videoSizeY = 480;

% define eye movements
% define rotations in deg
allPupilAzi = randn(1,nFrames)*15; % in degrees
allPupilEle = randn(1,nFrames)*15; % in degrees

%% construct scene geometry
% we will construct a scene where the CoR sits right at the center, and the
% pupil rotates on a 150 pixel radius.

% define the eye sphere radius
eyeRadius =  150;
pupilRadius = 30;

% create scene plane
sceneDistance = 1200; % orthogonal distance from rotation arm
imagePlaneEyeCenterX= videoSizeX/2; % max x size of scene (for plotting purposes)
imagePlaneEyeCenterY = videoSizeY/2; % max y size of scene

% Set up the sceneGeometry
projectionModel = 'pseudoPerspective';
eyeCenter = [imagePlaneEyeCenterX, imagePlaneEyeCenterY, eyeRadius+sceneDistance];
sceneGeometry.eyeCenter.X = imagePlaneEyeCenterX;
sceneGeometry.eyeCenter.Y = imagePlaneEyeCenterY;
sceneGeometry.eyeCenter.Z = eyeRadius+sceneDistance;
sceneGeometry.eyeRadius = eyeRadius;
sceneGeometry.meta.projectionModel = projectionModel;

% Set up anonymous function for nonlinear constraint
nonlinconst = @(transparentEllipseParams) constrainEllipseBySceneGeometry(...
    transparentEllipseParams, ...
    sceneGeometry, 1);

fprintf('Generating ellipses on the image plane...\n');

h = figure('visible','off');
emptyFrame = ones(videoSizeY, videoSizeX);

% for each frame, create the projected ellipse and test the consistency of
% the returned theta values
for ii = 1:nFrames
    
    % Produce the forward projection ellipse
    forwardProjectEllipseParams = pupilProjection_fwd(allPupilAzi(ii), allPupilEle(ii), pi*pupilRadius.^2, eyeCenter, eyeRadius, projectionModel);
    
    % make the plot
    imshow(emptyFrame);
    hold on
    pFitImplicit = ellipse_ex2im(ellipse_transparent2ex(forwardProjectEllipseParams));
    fh=@(x,y) pFitImplicit(1).*x.^2 +pFitImplicit(2).*x.*y +pFitImplicit(3).*y.^2 +pFitImplicit(4).*x +pFitImplicit(5).*y +pFitImplicit(6);
    fimplicit(fh,[1, videoSizeX, 1, videoSizeY],'Color', 'b','LineWidth',1);
    
    % get the perimeter
    tmpFrame = getframe(gca);
    [Yp, Xp] = ind2sub(size(squeeze(tmpFrame.cdata(:,:,1))),find(squeeze(tmpFrame.cdata(:,:,1))==0));
    perimeter.data{ii}.Yp = Yp;
    perimeter.data{ii}.Xp = Xp;

end

close(h);

fprintf('Checking theta consistency...\n');

% now loop through the frames and check out theta consistency
for ii = 1:nFrames

    Yp = perimeter.data{ii}.Yp;
    Xp = perimeter.data{ii}.Xp;

    % Obtain the forward projection ellipse params
    forwardProjectEllipseParams = pupilProjection_fwd(allPupilAzi(ii), allPupilEle(ii), pi*pupilRadius.^2, eyeCenter, eyeRadius, projectionModel);

    % Perform an unconstrained ellipse fit
    [unconstrainedEllipseParams, unconstrainedRMSE] = ...
        constrainedEllipseFit(Xp, Yp, ...
        [], ...
        [], ...
        []);
   
    % Perform a constrained ellipse fit
    [constrainedEllipseParams, constrainedRMSE] = ...
        constrainedEllipseFit(Xp, Yp, ...
        [0, 0, 300, 0, 0], ...
        [videoSizeX,videoSizeY,20000,1.0, pi], ...
        nonlinconst);
    
    % Assemble a vector of the theta results
    thetaVals = [forwardProjectEllipseParams(5) constrainedEllipseParams(5) ];
    
    % Obtain a list of pair-wise differences
    thetaDiscrepancies = abs(bsxfun(@minus,thetaVals,thetaVals'));
    
    % Handle the polar nature of theta. That is, a difference between two
    % thetas of 3.1 rads is in fact a polar difference of 0.04 rads
    thetaDiscrepancies=bsxfun(@min,thetaDiscrepancies,abs(thetaDiscrepancies-pi));
    
    % If the ellipse eccentricity is above threshold, test if there are
    % discrepancies above tolerance in the calculated theta value
    if unconstrainedEllipseParams(4) > eccentricityThresh
        if max(max(thetaDiscrepancies)) > thetaDifferenceTolerance
            
            % Calculate the RMSE for the forward projection ellipse
            [~, forwardProjectRMSE] = ...
                constrainedEllipseFit(Xp, Yp, ...
                forwardProjectEllipseParams, ...
                forwardProjectEllipseParams, ...
                []);

            % Report the info for the three fit approaches
            fprintf('unconstrained theta: %0.3f, RMSE: %0.3f, constraint: %0.3f \n', ...
                unconstrainedEllipseParams(5), ...
                unconstrainedRMSE, ...
                nonlinconst(unconstrainedEllipseParams));
            fprintf('forward projection theta: %0.3f, RMSE: %0.3f, constraint: %0.3f \n', ...
                forwardProjectEllipseParams(5), ...
                forwardProjectRMSE, ...
                nonlinconst(forwardProjectEllipseParams));
            fprintf('constrained theta: %0.3f, RMSE: %0.3f, constraint: %0.3f \n', ...
                constrainedEllipseParams(5), ...
                constrainedRMSE, ...
                nonlinconst(constrainedEllipseParams));
            
            % Make a figure illustrating the fits
            figure(1); clf;
            plot(Xp, Yp, '.b');
            xlim([0 videoSizeX]);
            ylim([0 videoSizeY]);
            hold on
            % Add the unconstrained fit
            pFitImplicit = ellipse_ex2im(ellipse_transparent2ex(unconstrainedEllipseParams));
            fh=@(x,y) pFitImplicit(1).*x.^2 +pFitImplicit(2).*x.*y +pFitImplicit(3).*y.^2 +pFitImplicit(4).*x +pFitImplicit(5).*y +pFitImplicit(6);
            fimplicit(fh,[1, videoSizeX, 1, videoSizeY],'Color', 'b','LineWidth',1);
            % Add the forward projection fit
            pFitImplicit = ellipse_ex2im(ellipse_transparent2ex(forwardProjectEllipseParams));
            fh=@(x,y) pFitImplicit(1).*x.^2 +pFitImplicit(2).*x.*y +pFitImplicit(3).*y.^2 +pFitImplicit(4).*x +pFitImplicit(5).*y +pFitImplicit(6);
            fimplicit(fh,[1, videoSizeX, 1, videoSizeY],'Color', 'g','LineWidth',1);
            % Add the constrained fit
            pFitImplicit = ellipse_ex2im(ellipse_transparent2ex(constrainedEllipseParams));
            fh=@(x,y) pFitImplicit(1).*x.^2 +pFitImplicit(2).*x.*y +pFitImplicit(3).*y.^2 +pFitImplicit(4).*x +pFitImplicit(5).*y +pFitImplicit(6);
            fimplicit(fh,[1, videoSizeX, 1, videoSizeY],'Color', 'r','LineWidth',1);
            legend('points','unconstrained','forward','constrained');
            
            pause
        end
    end

end % loop over frames

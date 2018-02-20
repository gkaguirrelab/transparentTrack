function [pupilEllipseOnImagePlane, imagePoints, sceneWorldPoints, eyeWorldPoints, pointLabels, nodalPointIntersectError] = pupilProjection_fwd(eyePose, sceneGeometry, rayTraceFuncs, varargin)
% Project the pupil circle to an ellipse on the image plane
%
% Syntax:
%  [pupilEllipseOnImagePlane, imagePoints, sceneWorldPoints, eyeWorldPoints, pointLabels] = pupilProjection_fwd(eyePoses, sceneGeometry, rayTraceFuncs)
%
% Description:
%   Given the sceneGeometry--and optionaly ray tracing functions through
%   the cornea--this routine simulates a circular pupil on a spherical eye
%   and then measures the parameters of the ellipse (in transparent format)
%   of the projection of the pupil to the image plane.
%
%   The forward model is a perspective projection of an anatomically
%   accurate eye, with points positioned behind the cornea subject
%   to refractive displacement. The projection incorporates the intrinsic
%   properties of the camera, including any radial lens distortion.
%
% Notes:
%   Rotations - Eye rotation is given as azimuth, elevation, and torsion in
%   degrees. These values correspond to degrees of rotation of the eye
%   relative to a head-fixed (extrinsic) coordinate frame. Note that this
%   is different from an eye-fixed (intrinsic) coordinate frame (such as
%   the Fick coordinate sysem). Azimuth, Elevation of [0,0] corresponds
%   to the position of the eye when a line that connects the center of
%   rotation of the eye with the center of the pupil is normal to the image
%   plane. Positive rotations correspond to rightward, upward, translation
%   of the pupil center in the image. Torsion of zero corresponds to the
%   torsion of the eye when it is in primary position.
%
%   Units - Eye rotations are in units of degrees. However, the units of
%   theta in the transparent ellipse parameters are radians. This is in
%   part to help us keep the two units separate conceptually.
%
% Inputs:
%   eyePose               - A 1x4 vector provides values for [eyeAzimuth,
%                           eyeElevation, eyeTorsion, pupilRadius].
%                           Azimuth, elevation, and torsion are in units of
%                           head-centered (extrinsic) degrees, and pupil
%                           radius in mm.
%   sceneGeometry         - A structure; described in createSceneGeometry
%   rayTraceFuncs         - A structure; described in assembleRayTraceFuncs
%
% Optional key/value pairs:
%  'fullEyeModelFlag'     - Logical. Determines if the full posterior and
%                           anterior chamber eye model will be created.
%  'nPupilPerimPoints'    - The number of points that are distributed
%                           around the pupil circle. A minimum of 5 is
%                           required to uniquely specify the image ellipse.
%  'nIrisPerimPoints'     - The number of points that are distributed
%                           around the iris circle. A minimum of 5 is
%                           required to uniquely specify the image ellipse.
%  'posteriorChamberEllipsoidPoints' - The number of points that are on
%                           each latitude line of the posterior chamber
%                           ellipsoid. About 30 makes a nice image.
%  'anteriorChamberEllipsoidPoints' - The number of points that are on
%                           each longitude line of the anterior chamber
%                           ellipsoid. About 30 makes a nice image.
%  'calcNodalIntersectError' - Logical. Controls if the nodal point
%                           intersection error is calculated. Set to false
%                           by default as this adds time to the forward
%                           model.
%
% Outputs:
%   pupilEllipseOnImagePlane - A 1x5 vector with the parameters of the
%                           pupil ellipse on the image plane cast in
%                           transparent form.
%   imagePoints           - An nx2 matrix that specifies the x, y location
%                           on the image plane for each of the eyeWorld
%                           points.
%   sceneWorldPoints      - An nx3 matrix of the coordinates of the
%                           points of the eye model in the sceneWorld
%                           coordinate frame. If fullEyeModel is set to
%                           false, then there will only be points for the
%                           pupil perimeter. If fullEyeModel is true, then
%                           the entire model of ~1000 points will be
%                           returned.
%   eyeWorldPoints        - An nx3 matrix of the coordinates of the
%                           points of the eye model in the eyeWorld
%                           coordinate frame.
%   pointsLabels          - An nx1 cell array that identifies each of the
%                           points, from the set {'pupilCenter',
%                           'irisCenter', 'rotationCenter',
%                           'posteriorChamber', 'irisPerimeter',
%                           'pupilPerimeter',
%                           'anteriorChamber','cornealApex'}.
%   nodalPointIntersectError - A nx1 vector that contains the distance (in
%                           mm) between the nodal point of the camera and
%                           the intersection of a ray on the camera plane.
%                           The value is nan for points not subject to
%                           refraction by the cornea. All values will be
%                           nan if the key value calcNodalIntersectError is
%                           set to false.
%
% Examples:
%{
    %% Obtain the parameters of the pupil ellipse in the image
    % Obtain a default sceneGeometry structure
    sceneGeometry=createSceneGeometry();
    % Define the ray tracing functions
    rayTraceFuncs = assembleRayTraceFuncs(sceneGeometry);
    % Define in eyePoses the azimuth, elevation, torsion, and pupil radius
    eyePose = [-10 5 0 3];
    % Obtain the pupil ellipse parameters in transparent format
    pupilEllipseOnImagePlane = pupilProjection_fwd(eyePose,sceneGeometry,rayTraceFuncs);
%}
%{
    %% Display a 2D image of a slightly myopic left eye
    % Obtain a default sceneGeometry structure
    sceneGeometry=createSceneGeometry('eyeLaterality','left','spectacleRefractionDiopters',-2);
    % Define the ray tracing functions
    rayTraceFuncs = assembleRayTraceFuncs(sceneGeometry);
    % Define an eyePose with azimuth, elevation, torsion, and pupil radius
    eyePose = [-10 5 0 3];
    % Perform the projection and request the full eye model
    [~, imagePoints, ~, ~, pointLabels] = pupilProjection_fwd(eyePose,sceneGeometry,rayTraceFuncs,'fullEyeModelFlag',true);
    % Define some settings for display
    eyePartLabels = {'rotationCenter', 'posteriorChamber' 'irisPerimeter' 'pupilPerimeter' 'anteriorChamber' 'cornealApex'};
    plotColors = {'+r' '.w' '.b' '*g' '.y' '*y'};
    blankFrame = zeros(480,640)+0.5;
    % Prepare a figure
    figure
    imshow(blankFrame, 'Border', 'tight');
    hold on
    axis off
    axis equal
    xlim([0 640]);
    ylim([0 480]);
    % Plot each anatomical component
    for pp = 1:length(eyePartLabels)
    	idx = strcmp(pointLabels,eyePartLabels{pp});
        plot(imagePoints(idx,1), imagePoints(idx,2), plotColors{pp})
    end
%}
%{
    %% Display a 3D plot of a right eye
    % Obtain a default sceneGeometry structure
    sceneGeometry=createSceneGeometry();
    % Define the ray tracing functions
    rayTraceFuncs = assembleRayTraceFuncs(sceneGeometry);
    % Define an eyePose with azimuth, elevation, torsion, and pupil radius
    eyePose = [0 0 0 3];
    % Perform the projection and request the full eye model
    [~, ~, ~, eyeWorldPoints, pointLabels] = pupilProjection_fwd(eyePose,sceneGeometry,rayTraceFuncs,'fullEyeModelFlag',true);
    % Define some settings for display
    eyePartLabels = {'rotationCenter', 'posteriorChamber' 'irisPerimeter' 'pupilPerimeter' 'anteriorChamber' 'cornealApex'};
    plotColors = {'+r' '.k' '*b' '*g' '.y' '*y'};
    % Prepare a figure
    figure
    % Plot each anatomical component
    for pp = 1:length(eyePartLabels)
    	idx = strcmp(pointLabels,eyePartLabels{pp});
        plot3(eyeWorldPoints(idx,1), eyeWorldPoints(idx,2), eyeWorldPoints(idx,3), plotColors{pp})
        hold on
    end
    hold off
    axis equal
%}


%% input parser
p = inputParser; p.KeepUnmatched = true;

% Required
p.addRequired('eyePose',@isnumeric);
p.addRequired('sceneGeometry',@(x)(isempty(x) || isstruct(x)));
p.addRequired('rayTraceFuncs',@(x)(isempty(x) || isstruct(x)));

p.addParameter('fullEyeModelFlag',false,@islogical);
p.addParameter('nPupilPerimPoints',5,@(x)(isnumeric(x) && x>=4));
p.addParameter('nIrisPerimPoints',5,@(x)(isnumeric(x) && x>=4));
p.addParameter('posteriorChamberEllipsoidPoints',30,@isnumeric);
p.addParameter('anteriorChamberEllipsoidPoints',30,@isnumeric);
p.addParameter('calcNodalIntersectError',false,@islogical);

% parse
p.parse(eyePose, sceneGeometry, rayTraceFuncs, varargin{:})


%% Check the input

if isempty(sceneGeometry)
    % No sceneGeometry was provided. Use the default settings
    sceneGeometry = createSceneGeometry();
end


%% Prepare variables
% Separate the eyePoses into individual variables
eyeAzimuth = eyePose(1);
eyeElevation = eyePose(2);
eyeTorsion = eyePose(3);
pupilRadius = eyePose(4);
nPupilPerimPoints = p.Results.nPupilPerimPoints;


%% Define an eye in eyeWorld coordinates
% This coordinate frame is in mm units and has the dimensions (p1,p2,p3).
% The diagram is of a cartoon pupil, being viewed directly from the front.
%
% Coordinate [0,0,0] corresponds to the apex (front surface) of the cornea,
% with the model eye having the property of the optical and pupil axes of
% the eye being aligned. The first dimension is depth, and has a negative
% value towards the back of the eye.
%
%                 |
%     ^         __|__
%  -  |        /     \
% p3  -  -----(   +   )-----
%  +  |        \_____/
%     v           |
%                 |
%
%           - <--p2--> +
%
% For the right eye, negative values on the p2 dimension are more temporal,
% and positive values are more nasal. Positive values of p3 are downward,
% and negative values are upward

% Define points around the pupil circle
perimeterPointAngles = 0:2*pi/nPupilPerimPoints:2*pi-(2*pi/nPupilPerimPoints);
eyeWorldPoints(1:nPupilPerimPoints,3) = ...
    sin(perimeterPointAngles)*pupilRadius + sceneGeometry.eye.pupilCenter(3);
eyeWorldPoints(1:nPupilPerimPoints,2) = ...
    cos(perimeterPointAngles)*pupilRadius + sceneGeometry.eye.pupilCenter(2);
eyeWorldPoints(1:nPupilPerimPoints,1) = ...
    0 + sceneGeometry.eye.pupilCenter(1);

% Create labels for the pupilPerimeter points
tmpLabels = cell(nPupilPerimPoints, 1);
tmpLabels(:) = {'pupilPerimeter'};
pointLabels = tmpLabels;

% If the fullEyeModel flag is set, then we will create a model of the
% posterior and anterior chambers of the eye.
if p.Results.fullEyeModelFlag
    
    % Add points for the center of the pupil, iris, and rotation
    eyeWorldPoints = [eyeWorldPoints; sceneGeometry.eye.pupilCenter];
    pointLabels = [pointLabels; 'pupilCenter'];
    eyeWorldPoints = [eyeWorldPoints; sceneGeometry.eye.irisCenter];
    pointLabels = [pointLabels; 'irisCenter'];
    eyeWorldPoints = [eyeWorldPoints; sceneGeometry.eye.rotationCenter];
    pointLabels = [pointLabels; 'rotationCenter'];
    
    % Create the posterior chamber ellipsoid. We switch dimensions here so
    % that the ellipsoid points have their poles at corneal apex and
    % posterior apex of the eye
    [p3tmp, p2tmp, p1tmp] = ellipsoid( ...
        sceneGeometry.eye.posteriorChamberCenter(3), ...
        sceneGeometry.eye.posteriorChamberCenter(2), ...
        sceneGeometry.eye.posteriorChamberCenter(1), ...
        sceneGeometry.eye.posteriorChamberRadii(3), ...
        sceneGeometry.eye.posteriorChamberRadii(2), ...
        sceneGeometry.eye.posteriorChamberRadii(1), ...
        p.Results.posteriorChamberEllipsoidPoints);
    % Convert the surface matrices to a vector of points and switch the
    % axes back
    ansTmp = surf2patch(p1tmp, p2tmp, p3tmp);
    posteriorChamberPoints=ansTmp.vertices;
    
    % Retain those points that are anterior to the center of the posterior
    % chamber and are posterior to the iris plane
    retainIdx = logical(...
        (posteriorChamberPoints(:,1) > sceneGeometry.eye.posteriorChamberCenter(1)) .* ...
        (posteriorChamberPoints(:,1) < sceneGeometry.eye.irisCenter(1)) ...
        );
    if all(~retainIdx)
        error('The iris center is behind the center of the posterior chamber');
    end
    posteriorChamberPoints = posteriorChamberPoints(retainIdx,:);
    
    % Add the points and labels
    eyeWorldPoints = [eyeWorldPoints; posteriorChamberPoints];
    tmpLabels = cell(size(posteriorChamberPoints,1), 1);
    tmpLabels(:) = {'posteriorChamber'};
    pointLabels = [pointLabels; tmpLabels];
    
    % Define points around the perimeter of the iris
    nIrisPerimPoints = p.Results.nIrisPerimPoints;
    perimeterPointAngles = 0:2*pi/nIrisPerimPoints:2*pi-(2*pi/nIrisPerimPoints);
    irisPoints(1:nIrisPerimPoints,3) = ...
        sin(perimeterPointAngles)*sceneGeometry.eye.irisRadius + sceneGeometry.eye.irisCenter(3);
    irisPoints(1:nIrisPerimPoints,2) = ...
        cos(perimeterPointAngles)*sceneGeometry.eye.irisRadius + sceneGeometry.eye.irisCenter(2);
    irisPoints(1:nIrisPerimPoints,1) = ...
        0 + sceneGeometry.eye.irisCenter(1);
    
    % Add the points and labels
    eyeWorldPoints = [eyeWorldPoints; irisPoints];
    tmpLabels = cell(size(irisPoints,1), 1);
    tmpLabels(:) = {'irisPerimeter'};
    pointLabels = [pointLabels; tmpLabels];
    
    % Create the anterior chamber ellipsoid.
    [p1tmp, p2tmp, p3tmp] = ellipsoid( ...
        sceneGeometry.eye.corneaFrontSurfaceCenter(1), ...
        sceneGeometry.eye.corneaFrontSurfaceCenter(2), ...
        sceneGeometry.eye.corneaFrontSurfaceCenter(3), ...
        sceneGeometry.eye.corneaFrontSurfaceRadius, ...
        sceneGeometry.eye.corneaFrontSurfaceRadius, ...
        sceneGeometry.eye.corneaFrontSurfaceRadius, ...
        p.Results.anteriorChamberEllipsoidPoints);
    % Convert the surface matrices to a vector of points and switch the
    % axes back
    ansTmp = surf2patch(p1tmp, p2tmp, p3tmp);
    anteriorChamberPoints=ansTmp.vertices;
    
    % Retain those points that are anterior to the iris plane
    retainIdx = logical(...
        (anteriorChamberPoints(:,1) > sceneGeometry.eye.irisCenter(1)));
    if all(~retainIdx)
        error('The pupil plane is set in front of the corneal apea');
    end
    anteriorChamberPoints = anteriorChamberPoints(retainIdx,:);
    
    % Add the points and labels
    eyeWorldPoints = [eyeWorldPoints; anteriorChamberPoints];
    tmpLabels = cell(size(anteriorChamberPoints,1), 1);
    tmpLabels(:) = {'anteriorChamber'};
    pointLabels = [pointLabels; tmpLabels];
    
    % Add a point for the corneal apex
    cornealApex=[0 0 0];
    eyeWorldPoints = [eyeWorldPoints; cornealApex];
    pointLabels = [pointLabels; 'cornealApex'];
end


%% Project the eyeWorld points to headWorld coordinates.
% This coordinate frame is in mm units and has the dimensions (h1,h2,h3).
% The diagram is of a cartoon eye, being viewed directly from the front.
%
%  h1 values negative --> towards the head, positive towards the camera
%
%         h2
%    0,0 ---->
%     |
%  h3 |
%     v
%
%               |
%             __|__
%            /  _  \
%    -------(  (_)  )-------  h2 (horizontal axis of the head)
%            \_____/          rotation about h2 causes pure vertical
%               |             eye movement
%               |
%
%               h3
%   (vertical axis of the head)
%  rotation about h3 causes pure
%     horizontal eye movement
%
%
%
% Position [0,-,-] indicates the front surface of the eye.
% Position [-,0,0] indicates the h2 / h3 position of the pupil axis of
% the eye when it is normal to the image plane.
%
% We will convert from this coordinate frame to that of the camera scene
% later.


%% Define the eye rotation matrix
% Assemble a rotation matrix from the head-fixed Euler angle rotations. In
% the head-centered world coordinate frame, positive azimuth, elevation and
% torsion values correspond to leftward, downward and clockwise (as seen
% from the perspective of the subject) eye movements
R3 = [cosd(eyeAzimuth) -sind(eyeAzimuth) 0; sind(eyeAzimuth) cosd(eyeAzimuth) 0; 0 0 1];
R2 = [cosd(eyeElevation) 0 sind(eyeElevation); 0 1 0; -sind(eyeElevation) 0 cosd(eyeElevation)];
R1 = [1 0 0; 0 cosd(eyeTorsion) -sind(eyeTorsion); 0 sind(eyeTorsion) cosd(eyeTorsion)];

% This order (1-2-3) corresponds to a head-fixed, extrinsic, rotation
% matrix. The reverse order (3-2-1) would be an eye-fixed, intrinsic
% rotation matrix and would corresponds to the "Fick coordinate" scheme.
eyeRotation = R1*R2*R3;


%% Obtain the virtual image for the eyeWorld points
% This steps accounts for the effect of corneal refraction upon the
% appearance of points from the iris and pupil
if ~isempty(rayTraceFuncs)
    % Identify the eyeWorldPoints that are subject to refraction by the cornea
    refractPointsIdx = find(strcmp(pointLabels,'pupilPerimeter')+...
        strcmp(pointLabels,'irisPerimeter')+...
        strcmp(pointLabels,'pupilCenter')+...
        strcmp(pointLabels,'irisCenter'));
    % Loop through the eyeWorldPoints that are to be refracted
    nodalPointIntersectError = nan(length(pointLabels),1);
    for ii=1:length(refractPointsIdx)
        % Grab this eyeWorld point
        eyeWorldPoint=eyeWorldPoints(refractPointsIdx(ii),:);
        % Define an error function which is the distance between the nodal
        % point of the camera and the point at which a ray impacts the
        % plane that contains the camera, with the ray departing from the
        % eyeWorld point at angle theta in the p1p2 plane.
        errorFunc = @(theta) rayTraceFuncs.cameraNodeDistanceError2D.p1p2(...
            sceneGeometry.extrinsicTranslationVector(1),...
            sceneGeometry.extrinsicTranslationVector(2),...
            sceneGeometry.extrinsicTranslationVector(3),...
            deg2rad(eyeAzimuth), deg2rad(eyeElevation), deg2rad(eyeTorsion),...
            eyeWorldPoint(1),eyeWorldPoint(2),eyeWorldPoint(3),...
            sceneGeometry.eye.rotationCenter(1),...
            theta);
        % Conduct an fminsearch to find the p1p2 theta that results in a
        % ray that strikes as close as possible to the camera nodal point.
        % Because the errorFunc returns nan for values very close to zero,
        % we initialize the search with a value slightly away (1e-4)
        theta_p1p2=fminsearch(errorFunc,1e-4);
        % Now repeat this process for a ray that varies in theta in the
        % p1p3 plane
        errorFunc = @(theta) rayTraceFuncs.cameraNodeDistanceError2D.p1p3(...
            sceneGeometry.extrinsicTranslationVector(1),...
            sceneGeometry.extrinsicTranslationVector(2),...
            sceneGeometry.extrinsicTranslationVector(3),...
            deg2rad(eyeAzimuth), deg2rad(eyeElevation), deg2rad(eyeTorsion),...
            eyeWorldPoint(1),eyeWorldPoint(2),eyeWorldPoint(3),...
            sceneGeometry.eye.rotationCenter(1),...
            theta);
        theta_p1p3=fminsearch(errorFunc,1e-4);
        % With both theta values calculated, now obtain the virtual image
        % ray arising from the pupil plane that reflects the corneal optics
        virtualImageRay = rayTraceFuncs.virtualImageRay(eyeWorldPoint(1), eyeWorldPoint(2), eyeWorldPoint(3), theta_p1p2, theta_p1p3);
        % Replace the original eyeWorld point with the virtual image
        % eyeWorld point
        eyeWorldPoints(refractPointsIdx(ii),:) = virtualImageRay(1,:);
        % The code below may be used to calculate the total error (in mm)
        % in both dimensions for intersecting the nodal point of the
        % camera. Error values on the order of 0.1 - 5 are found across
        % pupil points and for a range of eye rotations. By default, the
        % flag that controls this calculatuin is set to false, as the
        % computation is lengthy and is not otherwise used.
        if p.Results.calcNodalIntersectError
            nodalPointIntersectError(refractPointsIdx(ii)) = ...
                rayTraceFuncs.cameraNodeDistanceError3D(...
                sceneGeometry.extrinsicTranslationVector(1),...
                sceneGeometry.extrinsicTranslationVector(2),...
                sceneGeometry.extrinsicTranslationVector(3),...
                deg2rad(eyeAzimuth), deg2rad(eyeElevation), deg2rad(eyeTorsion),...
                eyeWorldPoint(1),eyeWorldPoint(2),eyeWorldPoint(3),...
                sceneGeometry.eye.rotationCenter(1),...
                theta_p1p2, theta_p1p3);
        end
    end
end

%% Apply the eye rotation
headWorldPoints = (eyeRotation*(eyeWorldPoints-sceneGeometry.eye.rotationCenter)')'+sceneGeometry.eye.rotationCenter;


%% Project the headWorld points to sceneWorld coordinates.
% This coordinate frame is in mm units and has the dimensions (X,Y,Z).
% The diagram is of a cartoon head (borrowed from Leszek Swirski), being
% viewed from above:
%
%   |
%   |    .-.
%   |   |   | <- Head
%   |   `^u^'
% Z |      :V <- Camera    (As seen from above)
%   |      :
%   |      :
%  \|/     o <- Target
%
%     ----------> X
%
% +X = right
% +Y = up
% +Z = front (towards the camera)
%
% The origin [0,0,0] corresponds to the front surface of the eye and the
% pupil center when the line that connects the center of rotation of the
% eye and the pupil center are normal to the image plane.

% Re-arrange the head world coordinate frame to transform to the scene
% world coordinate frame
sceneWorldPoints = headWorldPoints(:,[2 3 1]);


%% Project the sceneWorld points to the image plane
% This coordinate frame is in units of pixels, and has the dimensions
% [x, y]:
%
%      ^
%      |
%   y  |
%      |
%      +------->
% [0,0]    x
%
% With x being left/right and y being down/up
%


% Create the projectionMatrix
projectionMatrix = ...
    sceneGeometry.intrinsicCameraMatrix * ...
    [sceneGeometry.extrinsicRotationMatrix, ...
    sceneGeometry.extrinsicTranslationVector];

% What is our total number of points to project?
nEyeWorldPoints = size(eyeWorldPoints,1);

% Project the sceneWorld points to the image plane and scale. The
% sceneWorld points have a column of ones added support the multiplication
% with a combined rotation and translation matrix
tmpImagePoints=(projectionMatrix*[sceneWorldPoints, ones(nEyeWorldPoints,1)]')';
imagePointsPreDistortion=zeros(nEyeWorldPoints,2);
imagePointsPreDistortion(:,1) = ...
    tmpImagePoints(:,1)./tmpImagePoints(:,3);
imagePointsPreDistortion(:,2) = ...
    tmpImagePoints(:,2)./tmpImagePoints(:,3);


%% Apply radial lens distortion
% This step introduces "pincushion" (or "barrel") distortion produced by
% the lens. The x and y distortion equations are in the normalized image
% coordinates. Thus, the origin is at the sensor optical center (aka
% principal point), and the coordinates are in world units. To apply this
% distortion to our image coordinate points, we subtract the optical
% center, and then divide by fx and fy from the intrinsic matrix.
imagePointsNormalized = (imagePointsPreDistortion - [sceneGeometry.intrinsicCameraMatrix(1,3) sceneGeometry.intrinsicCameraMatrix(2,3)]) ./ ...
    [sceneGeometry.intrinsicCameraMatrix(1,1) sceneGeometry.intrinsicCameraMatrix(2,2)];

% Distortion is proportional to distance from the center of the center of
% projection on the camera sensor
radialPosition = sqrt(imagePointsNormalized(:,1).^2 + imagePointsNormalized(:,2).^2);

distortionVector =   1 + ...
    sceneGeometry.radialDistortionVector(1).*radialPosition.^2 + ...
    sceneGeometry.radialDistortionVector(2).*radialPosition.^4;

imagePointsNormalizedDistorted(:,1) = imagePointsNormalized(:,1).*distortionVector;
imagePointsNormalizedDistorted(:,2) = imagePointsNormalized(:,2).*distortionVector;

% Place the distorted points back into the imagePoints vector
imagePoints = (imagePointsNormalizedDistorted .* [sceneGeometry.intrinsicCameraMatrix(1,1) sceneGeometry.intrinsicCameraMatrix(2,2)]) +...
    [sceneGeometry.intrinsicCameraMatrix(1,3) sceneGeometry.intrinsicCameraMatrix(2,3)];


%% Fit an ellipse to the pupil points in the image plane
% Obtain the transparent ellipse params of the projection of the pupil
% circle on the image plane.
pupilPerimIdx = find(strcmp(pointLabels,'pupilPerimeter'));

% Before we try to fit the ellipse, make sure that the radius is not zero,
% that the image points are real (not imaginary points that did not pass
% through corneal optics), and that there are at least 5 perimeter points.
if eyePose(4)==0 || ~isreal(imagePoints(pupilPerimIdx,:)) || length(pupilPerimIdx)<5
    pupilEllipseOnImagePlane=nan(1,5);
else
    % We place the ellipse fit in a try-catch block, as the fit can
    % fail when the ellipse is so eccentric that it approaches a line
    try
        implicitEllipseParams = ellipsefit_direct( imagePoints(pupilPerimIdx,1), imagePoints(pupilPerimIdx,2));
        % Convert the ellipse from implicit to transparent form
        pupilEllipseOnImagePlane = ellipse_ex2transparent(ellipse_im2ex(implicitEllipseParams));
        % place theta within the range of 0 to pi
        if pupilEllipseOnImagePlane(5) < 0
            pupilEllipseOnImagePlane(5) = pupilEllipseOnImagePlane(5)+pi;
        end
    catch
        % In the event of an error, return nans for the ellipse
        pupilEllipseOnImagePlane = nan(1,length(pupilPerimIdx));
    end
end

end % pupilProjection_fwd




function opticalSystemOut = addSpectacle(opticalSystemIn, spectacleRefractionDiopters, varargin)
% Add a spectacle lens to a passed optical system
%
% Description:
%	This routine adds a meniscus (ophthalmologic) lens to an optical system
%	with the refractive power specified in the passed variable. The routine
%   assumes that the prior optical medium to encountering the surfaces of
%   the lens was air.
%
% Inputs:
%   opticalSystemIn       - An mx3 matrix, where m is the number of
%                           surfaces in the model, including the initial
%                           position of the ray. Each row contains the
%                           values [center, radius, refractiveIndex] that
%                           define a spherical lens.
%   spectacleRefractionDiopters - Scalar. Refractive power in units of 
%                           diopters. A negative value specifies a lens
%                           that would be worn by someone with myopia to
%                           correct their vision.
%
% Optional key/value pairs:
%  'vertexDistance'       - Scalar. Distance (in mm) between the corneal
%                           apexa and the back surface of the lens. Typical
%                           values are 12-14 mm/
%  'lensRefractiveIndex'  - Scalar. Refractive index of the lens material.
%                           Typical values are in the range of 1.5 - 1.7.
%  'nearPlanoCurvature'   - Scalar. This defines the curvature of the near-
%                           plano face of the lens.
%
% Outputs:
%   opticalSystemOut      - An (m+2)x3 matrix, correspinding to the
%                           opticalSystemIn with the addition of the
%                           spectacle lens
%
% Examples:
%{
    %% Example - Pupil through cornea and spectacle, plot range limits
    % A model of the passage of a point on the pupil perimeter through
    % the cornea and spectacle lens (units in mm)
    clear coords
    clear theta
    clear figureFlag
    %  Obtain the eye parameters from the modelEyeParameters() function
    eye = modelEyeParameters('spectacleRefractionDiopters',-2);
    pupilRadius = 2;
    theta = deg2rad(-45);
    coords = [eye.pupilCenter(1) pupilRadius];
    opticalSystem = [nan nan eye.aqueousRefractiveIndex; ...
                     eye.corneaBackSurfaceCenter(1) -eye.corneaBackSurfaceRadius eye.corneaRefractiveIndex; ...
                     eye.corneaFrontSurfaceCenter(1) -eye.corneaFrontSurfaceRadius 1.0];
    % Add a -2 diopter lens for the correction of myopia
    opticalSystem=addSpectacle(opticalSystem, 2);
    % Define FigureFlag as a structure, and set the new field to false so
    % that subsequent calls to the ray tracing routine will plot on the
    % same figure. Also, set the textLabels to false to reduce clutter
    figureFlag.show = true;
    figureFlag.new = true;
    figureFlag.surfaces = true;
    figureFlag.imageLines = true;
    figureFlag.rayLines = true;
    figureFlag.textLabels = true;
    figureFlag.zLim = [-20 20];
    figureFlag.hLim = [-25 25];
    outputRay = rayTraceCenteredSphericalSurfaces(coords, theta, opticalSystem, figureFlag)
%}

%% input parser
p = inputParser; p.KeepUnmatched = true;

% Required
p.addRequired('opticalSystemIn',@isnumeric);
p.addRequired('spectacleRefractionDiopters',@isnumeric);

% Optional
p.addParameter('vertexDistance',12,@isnumeric);
p.addParameter('lensRefractiveIndex',1.498,@isnumeric);
p.addParameter('nearPlanoCurvature',-16000,@isnumeric);

% parse
p.parse(opticalSystemIn, spectacleRefractionDiopters, varargin{:})

% Distribute the parameters into variables
vertexDistance = p.Results.vertexDistance;
lensRefractiveIndex = p.Results.lensRefractiveIndex;
nearPlanoCurvature = p.Results.nearPlanoCurvature;

% Copy the optical system from input to output
opticalSystemOut = opticalSystemIn;

% The lens equations do not perform properly for corrections of less that
% 0.25 diopters, and we don't bother trying to model so small a correction.
% In such a case, return the opticalSystem unaltered.
if abs(spectacleRefractionDiopters) < 0.25
    return
end

% Define a "meniscus" lens with two surfaces. Note that a ray emerging from
% the eye encounters two concave surfaces, so both will have a negative
% radius of curvature for rayTraceCenteredSphericalSurfaces()
if spectacleRefractionDiopters > 0
    % This is a plus lens for the correction of hyperopia. It has a
    % relatively flat back surface and a more curved front surface.
    % We first add the near-plano back surface to the optical system
    backCurvature = nearPlanoCurvature;
    backCenter = vertexDistance+backCurvature;
    opticalSystemOut(end+1,:)=[backCenter backCurvature lensRefractiveIndex];
    backDiopters = (1-lensRefractiveIndex)/(backCurvature/1000);
 
    % How many diopters of correction do we need from the front surface?
    frontDiopters = spectacleRefractionDiopters - backDiopters;
 
    % Calculate the radius of curvature of the front surface
    frontCurvature = ((1-lensRefractiveIndex)/frontDiopters)*1000;
    
    % Calculate the location of the center of curvature for the front lens.
    % To do so, we introduce lens thickness as a symbolic variable.
    syms thickness
    assume(thickness>0);
    frontCenter = vertexDistance + frontCurvature + thickness;
        
    % Define the equations for the height of each surface
    syms x
    assume(x>0);
    backPerimHeight = sqrt(backCurvature^2 - (x-backCenter)^2);
    frontPerimHeight = sqrt(frontCurvature^2 - (x-frontCenter)^2);
    
    % The lens will be thickest along the optical axis of the eye. We wish
    % the lens to have a non-negative thickness within the range of rays
    % we wish to model as coming from the eye. We consider a line with a
    % slope of 2 that originates at the corneal apex and find the height at
    % which it intersects the back surface.
    eqn= backPerimHeight == x*2;
    intersectHeight = eval(solve(eqn));

    % Now solve for the thickness
    eqn = subs(frontPerimHeight,x,intersectHeight) == subs(backPerimHeight,x,intersectHeight);
    thickness = min(eval(solve(eqn)));

    % clear the remaining symbolic params
    clear x
    
    % Store the lens front surface in the optical system
    opticalSystemOut(end+1,:)=[frontCurvature+vertexDistance+thickness frontCurvature 1.0];
else
    % This is a minus lens for the correction of myopia. It has a
    % relatively flat front surface and a more curved back surface. It will
    % be thinnest at the center of the lens on the optical axis. For this
    % model lens, we allow the thickness to become zero at this point
    % We first determine the properties of the the near-plano front surface
    frontCurvature = nearPlanoCurvature;
    frontDiopters = (1-lensRefractiveIndex)/(frontCurvature/1000);
 
    % How many diopters of correction do we need from the back surface?
    backDiopters = frontDiopters - spectacleRefractionDiopters;
 
    % Calculate the radius of curvature of the back surface
    backCurvature = ((1-lensRefractiveIndex)/backDiopters)*1000;
    
    % Calculate the locations of the center of curvature for the back and
    % front lens.
    backCenter = vertexDistance+backCurvature;
    frontCenter = vertexDistance+frontCurvature;
    
    % Add the surfaces to the optical system
    opticalSystemOut(end+1,:)=[backCenter backCurvature lensRefractiveIndex];
    opticalSystemOut(end+1,:)=[frontCenter frontCurvature 1.0];

end

end % function - addSpectacle
function opticalSystemOut = addContactLens(opticalSystemIn, lensRefractionDiopters, varargin)
% Add a contact lens to a passed optical system
%
% Description:
%	This routine adds a meniscus (ophthalmologic) contact lens to an
%	optical system with the refractive power specified in the passed
%	variable. Note that a ray emerging from the eye encounters two concave
%	surfaces for this lens, so both surfaces will have a negative radius of
%	curvature for rayTraceCenteredSphericalSurfaces().
%
% Inputs:
%   opticalSystemIn       - An mx3 matrix, where m is the number of
%                           surfaces in the model, including the initial
%                           position of the ray. Each row contains the
%                           values [center, radius, refractiveIndex] that
%                           define a spherical lens.
%   lensRefractionDiopters - Scalar. Refractive power in units of
%                           diopters. A negative value specifies a lens
%                           that would be worn by someone with myopia to
%                           correct their vision.
%
% Optional key/value pairs:
%  'lensRefractiveIndex'  - Scalar. Refractive index of the lens material.
%                           Contact lens material is ~1.4.
%
% Outputs:
%   opticalSystemOut      - An (m+2)x3 matrix, correspinding to the
%                           opticalSystemIn with the addition of the
%                           contact lens
%
% Examples:
%{
    %% Example - Pupil through cornea and contact lens, plot range limits
    % A model of the passage of a point on the pupil perimeter through
    % the cornea and a contact lens (units in mm)
    clear coords
    clear theta
    clear figureFlag
    %  Obtain the eye parameters from the modelEyeParameters() function
    eye = modelEyeParameters('spectacleRefractionDiopters',-4);
    pupilRadius = 2;
    theta = deg2rad(-45);
    coords = [eye.pupilCenter(1) pupilRadius];
    opticalSystem = [nan nan eye.aqueousRefractiveIndex; ...
                     eye.corneaBackSurfaceCenter(1) -eye.corneaBackSurfaceRadius eye.corneaRefractiveIndex; ...
                     eye.corneaFrontSurfaceCenter(1) -eye.corneaFrontSurfaceRadius 1.0];
    % Add a -2 diopter lens for the correction of myopia
    opticalSystem=addContactLens(opticalSystem, -4);
    % Define FigureFlag as a structure, and set the new field to false so
    % that subsequent calls to the ray tracing routine will plot on the
    % same figure. Also, set the textLabels to false to reduce clutter
    figureFlag.show = true;
    figureFlag.new = true;
    figureFlag.surfaces = true;
    figureFlag.imageLines = true;
    figureFlag.rayLines = true;
    figureFlag.textLabels = true;
    figureFlag.zLim = [-10 5];
    figureFlag.hLim = [-5 5];
    outputRay = rayTraceCenteredSphericalSurfaces(coords, theta, opticalSystem, figureFlag)
%}

%% input parser
p = inputParser; p.KeepUnmatched = true;

% Required
p.addRequired('opticalSystemIn',@isnumeric);
p.addRequired('lensRefractionDiopters',@isnumeric);

% Optional
p.addParameter('lensRefractiveIndex',1.43,@isnumeric);
p.addParameter('minimumLensThickness',0.05,@isnumeric);

% parse
p.parse(opticalSystemIn, lensRefractionDiopters, varargin{:})

% Distribute the parameters into variables
lensRefractiveIndex = p.Results.lensRefractiveIndex;

% Copy the optical system from input to output
opticalSystemOut = opticalSystemIn;

% The lens equations do not perform properly for corrections of less that
% 0.25 diopters, and we don't bother trying to model so small a correction.
% In such a case, return the opticalSystem unaltered.
if abs(lensRefractionDiopters) < 0.25
    return
end

% The passed optical system will have a ray that emerges into a medium with
% a specified index of refraction. Because the contact lens contacts the
% corneal surface, this index of refraction will be replaced with the index
% of refraction of the contact lens. We store the index of refraction of
% the ambient medium (which will typically be air and thus 1.0) to apply to
% the final exit ray.
priorRefractiveIndex = opticalSystemIn(end-1,3);
mediumRefractiveIndex = opticalSystemIn(end,3);
opticalSystemOut(end,3) = lensRefractiveIndex;

% Calculate the diopters of the corneal surface without a contact lens; our
% goal is to create a front surface of the contact lens that produces a
% refractive correction equal to:
%   cornealSurfaceDiopters + lensRefractionDiopters
cornealSurfaceDiopters = (mediumRefractiveIndex-priorRefractiveIndex)/(opticalSystemIn(end,2)/1000);

% The back surface of the contact lens has the same radius of curvature
% and position as the front surface of the cornea. We then calculate
% the refractive power of the back surface of the lens, which will be
% small given that the index of refraction of the lens will be similar
% to that of the cornea.
backCurvature = opticalSystemIn(end,2);
backCenter = opticalSystemIn(end,1);
opticalSystemOut(end+1,:)=[backCenter backCurvature lensRefractiveIndex];
backDiopters = (lensRefractiveIndex-priorRefractiveIndex)/(backCurvature/1000);


% We calculate here thickness and thus center of the front surface.
if lensRefractionDiopters > 0
    % This is a plus lens for the correction of hyperopia.
    
    % How many diopters of correction do we need from the front surface?
    frontDiopters = lensRefractionDiopters - backDiopters + cornealSurfaceDiopters;
    
    % Calculate the radius of curvature of the front surface
    frontCurvature = ((mediumRefractiveIndex-lensRefractiveIndex)/frontDiopters)*1000;

    % Calculate the location of the center of curvature for the front lens.
    frontCenter = frontCurvature + (frontCurvature-backCurvature);
    
    % Store the lens front surface in the optical system
    opticalSystemOut(end+1,:)=[frontCenter frontCurvature mediumRefractiveIndex];
else
    % This is a minus lens for the correction of myopia.
    % It will be thinnest at the center of the lens on the optical axis.
    % The parameter minimumLensThickness defines this minimum.
    
    % How many diopters of correction do we need from the front surface?
    frontDiopters = lensRefractionDiopters - backDiopters + cornealSurfaceDiopters;
    
    % Calculate the radius of curvature of the front surface
    frontCurvature = ((mediumRefractiveIndex-lensRefractiveIndex)/frontDiopters)*1000;
    
    % Calculate the location of the center of curvature for the front lens.
    frontCenter = frontCurvature + p.Results.minimumLensThickness;

    % Add the surfaces to the optical system
    opticalSystemOut(end+1,:)=[frontCenter frontCurvature mediumRefractiveIndex];
    
end

end % function - addContactLens
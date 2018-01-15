function eye = modelEyeParameters( spectacleRefractionDiopters )
% Return the parameters of a model eye
%
% Description:
%   This routine returns the parameters that define the model eye that is
%   used in the sceneGeometry routines. We primarily make use of the values
%   derived by Atchison:
%
%       Atchison, David A. "Optical models for human myopic eyes." Vision 
%       research 46.14 (2006): 2236-2250.
%
%       Atchison, David A., et al. "Shape of the retinal surface in 
%       emmetropia and myopia." IOVS 46.8 (2005): 2698-2707.
%
%   Achison uses the dimensions [x, y, z] corresponding to the width, 
%   height, and depth (axial length) of the model eye. The parameters
%   returned by this routine correspond to the eyeWorld coordinate space
%   used in pupilProjection_fwd. The coordinates space is relative to the
%   apex of the cornea, with the dimensions [depth, vertical, horizontal],
%   and negative values of depth towards the center of the eye.
%
% Inputs:
%   spectacleRefractionDiopters - Scalar, in units of diopters. The
%                           dimensions of the posterior chamber of the eye
%                           (and to a lesser extent the curvature of the
%                           cornea) change with the observed refractive
%                           error of the subject. This value is the
%                           spherical refractive correction for the
%                           subject. A negative number is the correction
%                           that would be used for a myopic person.
%
% Outputs:
%   eye                   - A structure with fields that contain the values
%                           for the model eye.

% These values taken from Achison 2006, Table 1. The center of the cornea
% circle for the back surface is positioned so that there is 0.55 between
% the front and back surface of the cornea at the apex.
eye.corneaFrontSurfaceRadius = 7.77 + 0.022*spectacleRefractionDiopters;
eye.corneaFrontSurfaceCenter = [-7.77 0 0];
eye.corneaBackSurfaceRadius = 6.4;
eye.corneaBackSurfaceCenter = [-7.22 0 0];

% We position the pupil plane at the depth of the anterior point of the
% lens
eye.pupilCenter = [-3.7 0 0];

% Define the properties of the iris
eye.irisRadius = 6.5;
eye.irisCenter = [-4.5 0 0];

% The posterior chamber of the eye is modeled as an ellipsoid. The
% dimensions of the ellipsoid are set by the values observed for the radius
% of curvature of the retina in Atchison 2005. Note the greater dependence
% of the first dimension (axial length) upon refractive error.
eye.posteriorChamberRadii = [...
    10.148 - 0.163 * spectacleRefractionDiopters ...
    11.365 - 0.090 * spectacleRefractionDiopters ...
    11.455 - 0.043 * spectacleRefractionDiopters ];

% We position the posterior chamber by placing the anterior apex of the
% ellipsoid at a fixed distance behind the corneal apex. The value for this
% fixed distance is derived from the Atchison model eye as follows: The
% total interior depth of the eye from Table 1 for an emmetrope is:
%
%   0.55 + 3.15 + 1.44 + 2.16 + 16.28 = 23.5800
%
% The axial length of the posterior chamber ellipsoid in an emmetrope is:
%
%   10.1418 * 2 = 20.2836
%
% Therefore, the apex of the posterior ellipsoid is 23.5800 - 20.2836 =
% 3.2964 mm behind the corneal apex.
eye.posteriorChamberCenter = ...
    [(-3.2964 - eye.posteriorChamberRadii(1)) 0 0];

% The eye center of rotation is on average 13.3 mm behind the corneal
% apex per Gunter K. vonNoorden, MD; Emilio C. Campos "Binocular Vision and
% Ocular Motility Theory and Management of Strabismus" American Orthoptic
% Journal 51.1 (2001): 161-162.
eye.centerOfRotation = [-13.3 0 0];

% Refractive index values from Atchison 2006.
eye.corneaRefractiveIndex = 1.376;
eye.aqueousRefractiveIndex = 1.3374;

% Meta data regarding the units of the model
eye.meta.spectacleRefractionDiopters = spectacleRefractionDiopters;
eye.meta.units = 'mm';
eye.meta.coordinates = 'eyeWorld';
eye.meta.dimensions = {'depth (axial)' 'vertical' 'horizontal'};

end


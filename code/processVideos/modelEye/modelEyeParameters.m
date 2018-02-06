function eye = modelEyeParameters( varargin )
% Return the parameters of a model eye
%
% Syntax:
%  eye = modelEyeParameters()
%
% Description:
%   This routine returns the parameters that define the model eye that is
%   used in the sceneGeometry routines. We make use of the values derived
%   by Atchison for human eyes:
%
%       Atchison, David A. "Optical models for human myopic eyes." Vision
%       research 46.14 (2006): 2236-2250.
%
%       Atchison, David A., et al. "Shape of the retinal surface in
%       emmetropia and myopia." IOVS 46.8 (2005): 2698-2707.
%
%   Atchison uses the dimensions [x, y, z] corresponding to the width,
%   height, and depth (axial length) of the model eye. The parameters
%   returned by this routine correspond to the eyeWorld coordinate space
%   used in pupilProjection_fwd, which is relative to the apex of the
%   cornea and has the dimensions [depth, horizontal, vertical]; negative
%   values of depth towards the center of the eye.
%
% Inputs:
%   none
%
% Optional key/value pairs:
%  'spectacleRefractionDiopters' - Scalar, in units of diopters. The
%                           dimensions of the posterior chamber of the eye
%                           (and to a lesser extent the curvature of the
%                           cornea) change with the observed refractive
%                           error of the subject. This value is the
%                           spherical refractive correction for the
%                           subject. A negative number is the correction
%                           that would be used for a myopic person.
%  'axialLength'          - When set, this fixes the axial length of the 
%                           eye to the passed value in millimeters. As the
%                           modeled anterior chamber depth is not variable,
%                           this change is enforced on the posterior
%                           chamber. The remaining dimensions of the
%                           posterior chamber are scaled to fit the
%                           proportions predicted by the Atchison model for
%                           the specified degree of ametropia.
%  'eyeLaterality'        - A text string that specifies which eye (left,
%                           right) to model. Allowed values (in any case)
%                           are {'left','right','L','R','OS','OD'}
%  'species'              - A text string that specifies the species to be
%                           modeled. Supported values (in any case) are
%                           {'human'}
%                       
%
% Outputs:
%   eye                   - A structure with fields that contain the values
%                           for the model eye.
%
% Examples:
%{
    % Default parameters, corresponding to an emmetropic, right, human eye
    eye = modelEyeParameters();
%}
%{
    % Parameters for an myopic (-3), left, human eye
    eye = modelEyeParameters('spectacleRefractionDiopters',-3,'eyeLaterality','left');
%}


%% input parser
p = inputParser; p.KeepUnmatched = true;

% Optional
p.addParameter('spectacleRefractionDiopters',0,@isnumeric);
p.addParameter('axialLength',[],@(x)(isempty(x) || isnumeric(x)));
p.addParameter('eyeLaterality','Right',@ischar);
p.addParameter('species','Human',@ischar);

% parse
p.parse(varargin{:})

% Switch parameters at the top level by species
switch p.Results.species
    case {'human','Human','HUMAN'}

        % These values taken from Atchison 2006, Table 1. The center of the
        % cornea circle for the back surface is positioned so that there is
        % 0.55 between the front and back surface of the cornea at the
        % apex.
        eye.corneaFrontSurfaceRadius = 7.77;
        eye.corneaFrontSurfaceCenter = [-7.77 0 0];
        eye.corneaBackSurfaceRadius = 6.4;
        eye.corneaBackSurfaceCenter = [-7.22 0 0];
        
        % We position the pupil plane at the depth of the anterior point of
        % the lens. We have calculated for the optical model of the eye the
        % decentration values that produce an entrance pupil width / height
        % ratio that best matches empirical measurements (see
        % TEST_entrancePupilShape). The values from that search are used
        % here for the p2 and p3 position of the center of the pupil. In
        % both eyes, the pupil center is shifted slightly nasally and
        % downward. This is in keeping with empirical observations that:
        %
        %   ...the typical entrance pupil is decentered
        %   approximately 0.15 mm nasally and 0.1 mm inferior to the
        %   geometric center of the visible iris circumference
        %
        % Bennett, Edward S., and Barry A. Weissman, eds. Clinical contact
        % lens practice. Lippincott Williams & Wilkins, 2005, p119
        %
        switch p.Results.eyeLaterality
            case {'right','RIGHT','Right','R','r','od','OD'}
                eye.pupilCenter = [-3.7 0.15 0.1];
            case {'left','LEFT','Left','L','l','os','OS'}
                eye.pupilCenter = [-3.7 -0.15 0.1];
            otherwise
                error('Please specify a valid eye laterality for the model eye');
        end
        
        % Define the properties of the iris. Geoff needs to do some work to
        % find and justify values here.
        eye.irisRadius = 5;
        eye.irisCenter = [-4.5 0 0];
        
        % The posterior chamber of the eye is modeled as an ellipsoid. The
        % dimensions of the ellipsoid are set by the values observed for
        % the radius of curvature of the retina in Atchison 2005. Note the
        % greater dependence of the first dimension (axial length) upon
        % refractive error.        
        eye.posteriorChamberRadii = [...
            10.148 - 0.163 * p.Results.spectacleRefractionDiopters ...
            11.365 - 0.090 * p.Results.spectacleRefractionDiopters ...
            11.455 - 0.043 * p.Results.spectacleRefractionDiopters ];
        
        if ~isempty(p.Results.axialLength)
            % If a specific axial length is passed (perhaps obtained by
            % measurement using the IOL Master apparatus), set the model
            % eye to have this length, and scale the other dimensions of
            % the posterior chamber to maintain the specified ametropia. We
            % adjust the axial length for the component of the anterior
            % chamber that contibutes to length, which is calculated below
            % to be 3.2964.
            scaleFactor = (p.Results.axialLength - 3.2964) / (eye.posteriorChamberRadii(1)*2);
            eye.posteriorChamberRadii = eye.posteriorChamberRadii .* scaleFactor;
        end
        
        % We position the posterior chamber by placing the anterior apex of
        % the ellipsoid at a fixed distance behind the corneal apex. The
        % value for this fixed distance is derived from the Atchison model
        % eye as follows: The total interior depth of the eye from Table 1
        % for an emmetrope is:
        %
        %   0.55 + 3.15 + 1.44 + 2.16 + 16.28 = 23.5800
        %
        % The axial length of the posterior chamber ellipsoid in an
        % emmetrope is:
        %
        %   10.1418 * 2 = 20.2836
        %
        % Therefore, the apex of the posterior ellipsoid is 23.5800 -
        % 20.2836 = 3.2964 mm behind the corneal apex.
        eye.posteriorChamberCenter = ...
            [(-3.2964 - eye.posteriorChamberRadii(1)) 0 0];
        
        % The eye center of rotation is on average 13.3 mm behind the
        % corneal apex per Gunter K. vonNoorden, MD; Emilio C. Campos
        % "Binocular Vision and Ocular Motility Theory and Management of
        % Strabismus" American Orthoptic Journal 51.1 (2001): 161-162.
        % Spectacle refraction adjusts the axial length of the eye. We
        % assume here that the center of rotation reflects this change in
        % length. Specifically, we assume that the increase in the radius
        % of curvature of the posterior chamber described by the Atchison
        % model results in an equivalent lengthening of the radius of
        % rotation. Support for this 1:1 relationship is found in the paper
        %
        %   Dick, Graham L., Bryan T. Smith, and Peter L. Spanos. 
        %   "Axial length and radius of rotation of the eye." 
        %   Clinical and Experimental Optometry 73.2 (1990): 43-50.
        %
        % Specifically, Figure 6 shows that, for each mm of increase in 
        % the axial length of an eye, the center of rotation tended to
        % increase by 0.5 mm. Thus, there is a 1:1 relationship of axial
        % radius and rotation length.
        eye.rotationCenter = [(-3.152 - eye.posteriorChamberRadii(1)) 0 0];
        
        % Refractive index values from Atchison 2006.
        eye.corneaRefractiveIndex = 1.376;
        eye.aqueousRefractiveIndex = 1.3374;
        
    otherwise
        error('Please specify a valid species for the eye model');
end

% Meta data regarding the units of the model
eye.meta.spectacleRefractionDiopters = p.Results.spectacleRefractionDiopters;
eye.meta.axialLength = p.Results.axialLength;
eye.meta.laterality = p.Results.eyeLaterality;
eye.meta.species = p.Results.species;
eye.meta.units = 'mm';
eye.meta.coordinates = 'eyeWorld';
eye.meta.dimensions = {'depth (axial)' 'horizontal' 'vertical'};

end % function


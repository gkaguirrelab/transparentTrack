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
%   used in pupilProjection_fwd, which is relative to the pupil axis, with
%   the apex of the cornea set as zero in depth. The space has the
%   dimensions [depth, horizontal, vertical]; negative values of depth
%   towards the center of the eye. The model assumes the optical and pupil
%   axis of the eye are algined.
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
%  'axialLength'          - Scalar. When set, this fixes the axial length 
%                           of the eye to the passed value in millimeters.
%                           As the modeled anterior chamber depth is not
%                           variable, this change is enforced on the
%                           posterior chamber. The remaining dimensions of
%                           the posterior chamber are scaled to fit the
%                           proportions predicted by the Atchison model for
%                           the specified degree of ametropia.
%  'kappaAngle'           - 1x2 matrix. This is the angle of the visual
%                           axis in degrees w.r.t. to pupil axis. The
%                           values are [azimuth, elevation]. An eyePose of:
%                               [-kappa(1), -kappa(2), 0, radius]
%                           aligns the visual axis of the eye with the
%                           optical axis of the camera
%  'eyeLaterality'        - A text string that specifies which eye (left,
%                           right) to model. Allowed values (in any case)
%                           are {'left','right','L','R','OS','OD'}
%  'species'              - A text string that specifies the species to be
%                           modeled. Supported values (in any case) are
%                           {'human'}
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
p.addParameter('kappaAngle',[],@(x)(isempty(x) || isnumeric(x)));
p.addParameter('eyeLaterality','Right',@ischar);
p.addParameter('species','Human',@ischar);

% parse
p.parse(varargin{:})

% Interpret the passed laterality
switch p.Results.eyeLaterality
    case {'right','RIGHT','Right','R','r','od','OD'}
        eyeLaterality = 'Right';
    case {'left','LEFT','Left','L','l','os','OS'}
        eyeLaterality = 'Left';
    otherwise
        error('Please specify a valid eye laterality for the model eye');
end

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
        % the lens. The coordinate space of the model eye is define with
        % respect to the center of the pupil, so the p2 and p3 values are
        % zero
        eye.pupilCenter = [-3.7 0 0];
        
        % Define the properties of the iris. Geoff needs to do some work to
        % find and justify values here.
        % https://www.clspectrum.com/issues/2002/april-2002/contact-lens-case-reports        
        eye.irisRadius = 5;
                
        % In both eyes, the iris center is shifted slightly temporally and
        % upward with respect to the pupil center:
        %
        %   ...the typical entrance pupil is decentered
        %   approximately 0.15 mm nasally and 0.1 mm inferior to the
        %   geometric center of the visible iris circumference
        %
        % Bennett, Edward S., and Barry A. Weissman, eds. Clinical contact
        % lens practice. Lippincott Williams & Wilkins, 2005, p119
        switch eyeLaterality
            case 'Right'
                eye.irisCenter = [-4.5 -0.15 -0.1];
            case 'Left'
                eye.irisCenter = [-4.5 0.15 -0.1];
        end
        
        % The posterior chamber of the eye is modeled as an ellipsoid. The
        % dimensions of the ellipsoid are set by the values observed for
        % the radius of curvature of the retina in Atchison 2005. Note the
        % greater dependence of the first (axial) dimension upon
        % refractive error.        
        eye.posteriorChamberRadii = [...
            10.148 - 0.163 * p.Results.spectacleRefractionDiopters ...
            11.365 - 0.090 * p.Results.spectacleRefractionDiopters ...
            11.455 - 0.043 * p.Results.spectacleRefractionDiopters ];
        
        % Our model holds the depth of the anterior chamber constant. To
        % position the posterior chamber, we need to know the distance
        % between the apex of the anterior chamber and the apex of the
        % posterior chamber. The value for this fixed distance is derived
        % from the Atchison 2006 model eye. The total interior depth of the
        % eye from Table 1 for an emmetrope is:
        %
        %   0.55 + 3.15 + 1.44 + 2.16 + 16.28 = 23.58
        %
        % The axial length of the posterior chamber ellipsoid in an
        % emmetrope is:
        %
        %   10.1418 * 2 = 20.2836
        %
        % Therefore, the apex of the posterior ellipsoid is:
        %   23.5800 - 20.2836 = 3.2964 mm
        % behind the corneal apex.
        posteriorChamberApexDepth = 3.2964;
        
        % Compute and store axial length
        if isempty(p.Results.axialLength)
            eye.axialLength = posteriorChamberApexDepth + eye.posteriorChamberRadii(1)*2;
        else
            % If a specific axial length was passed (perhaps obtained by
            % measurement using the IOL Master apparatus), set the model
            % eye to have this length, and scale the other dimensions of
            % the posterior chamber to maintain the specified ametropia. We
            % adjust the axial length for the component of the anterior
            % chamber that contibutes to length, which is calculated below
            % to be 3.2964.
            scaleFactor = (p.Results.axialLength - posteriorChamberApexDepth) / (eye.posteriorChamberRadii(1)*2);
            eye.posteriorChamberRadii = eye.posteriorChamberRadii .* scaleFactor;
            eye.axialLength = p.Results.axialLength;
        end
        
        % Set the depth of the center of the posterior chamber
        eye.posteriorChamberCenter = ...
            [(-posteriorChamberApexDepth - eye.posteriorChamberRadii(1)) 0 0];
        
        % The eye center of rotation in emmetropes is on average 13.3 mm
        % behind the corneal apex per Gunter K. vonNoorden, MD; Emilio C.
        % Campos "Binocular Vision and Ocular Motility Theory and
        % Management of Strabismus" American Orthoptic Journal 51.1 (2001):
        % 161-162. Spectacle refraction adjusts the axial length of the
        % eye. We assume here that the center of rotation reflects this
        % change in length. Specifically, we assume that the increase in
        % the radius of curvature of the posterior chamber described by the
        % Atchison model produces an equivalent lengthening of the radius
        % of rotation. Support for this 1:1 relationship is found in the
        % paper:
        %
        %   Dick, Graham L., Bryan T. Smith, and Peter L. Spanos. 
        %   "Axial length and radius of rotation of the eye." 
        %   Clinical and Experimental Optometry 73.2 (1990): 43-50.
        %
        % Specifically, Figure 6 shows that, for each mm of increase in 
        % the axial length of an eye, the center of rotation tended to
        % increase by 0.5 mm. Thus, there is a 1:1 relationship of axial
        % radius and rotation length.
        % 
        % In an emmetropic eye, the distance from the corneal apex to the
        % center of the posterior chamber is
        %
        %   posteriorChamberApexDepth + posteriorChamberRadii(1), or
        %   3.2964 + 10.148 = 13.4444
        %
        % This implies that the center of rotation of an emmetropic eye
        % lies (13.4444 - 13.3 = 0.1444 mm) closer to the corneal surface
        % than the position of the center of the posterior chamber.
        eye.rotationCenter = [eye.posteriorChamberCenter(1)+0.1444 0 0];
        
        % Refractive index values from Atchison 2006.
        eye.corneaRefractiveIndex = 1.376;
        eye.aqueousRefractiveIndex = 1.3374;
        
        % We now calculate kappa, which is the angle (in degrees) between
        % the pupil and visual axes of the eye. The visual axis is
        % displaced nasally and superiorly within the visual field relative
        % to the pupil axis. Horizontal kappa is usually defined with
        % positive values being more nasal. We adopt a different convention
        % in which kappa is defined in head-fixed coordinates. Thus,
        % positive values for the right eye, and negative values for the
        % left eye, are more nasal. Positive values for vertical kappa are
        % upward.
        %
        % A source for an estimate of kappa comes from Mathur 2013:
        %
        %	Mathur, Ankit, Julia Gehrmann, and David A. Atchison. "Pupil shape
        %	as viewed along the horizontal visual field." Journal of vision
        %	13.6 (2013): 3-3.
        %
        % They measured the shape of the entrance pupil as a function of
        % viewing angle relative to the fixation point of the eye. Their
        % data is well fit by a kappa of [5, -2] degrees (see
        % TEST_entrancePupilShape.m). 
        %
        % Measured kappa has been found to depend upon axial length:
        %
        %   Tabernero, Juan, et al. "Mechanism of compensation of
        %   aberrations in the human eye." JOSA A 24.10 (2007): 3274-3283.
        %
        % Tabernero 2007 reports a mean horizontal kappa of 5 degrees in
        % emmetropes, and their Equation 6 expresses kappa (technically
        % alpha, the angle w.r.t. the optical axis) as a function of axial
        % length. Their formula assumes an emmetropic model eye of 24 mm,
        % while the model eye used here has an emmetropic axial length of
        % 23.592. The equation implemented below is adjusted so that an
        % emmetropic eye of 23.5924 mm has a horizontal (nasal directed)
        % kappa of 5 degrees and a vertical (inferiorly directed) kappa of
        % -2 degrees.
        %
        % While a horizontal kappa of ~5 degrees is a consistent finding,
        % measurements of vertical kappa differ: 
        %
        %   Hashemi, Hassan, et al. "Distribution of angle kappa
        %   measurements with Orbscan II in a population-based survey."
        %   Journal of Refractive Surgery 26.12 (2010): 966-971.
        %
        %   Gharaee, Hamid, et al. "Angle kappa measurements: normal values
        %   in healthy iranian population obtained with the Orbscan II."
        %   Iranian Red Crescent Medical Journal 17.1 (2015).
        %
        % We note that there is evidence that the vertical kappa value can
        % vary based upon the subject being in a sittng or supine position.
        % Until better evidene is available, we adopt a vertical kappa of
        % -2 degrees for the emmetropic model eye.

        if isempty(p.Results.kappaAngle)
            switch eyeLaterality
                case 'Right'
                    eye.kappaAngle(1) = atand((15.0924/(eye.axialLength-8.5000))*tand(5));
                case 'Left'
                    eye.kappaAngle(1) = -atand((15.0924/(eye.axialLength-8.5000))*tand(5));
            end
            eye.kappaAngle(2) = atand((15.0924/(eye.axialLength-8.5000))*tand(2));
        else
            eye.kappaAngle = p.Results.kappaAngle;
        end
        
    otherwise
        error('Please specify a valid species for the eye model');
end

% Meta data regarding the units of the model
eye.meta.spectacleRefractionDiopters = p.Results.spectacleRefractionDiopters;
eye.meta.axialLength = p.Results.axialLength;
eye.meta.laterality = eyeLaterality;
eye.meta.species = p.Results.species;
eye.meta.units = 'mm';
eye.meta.coordinates = 'eyeWorld';
eye.meta.dimensions = {'depth (axial)' 'horizontal' 'vertical'};
eye.meta.kappa = 'Degrees angle of visual axis w.r.t. pupil axis.';

end % function


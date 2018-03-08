function eye = modelEyeParameters( varargin )
% Return the parameters of a model eye
%
% Syntax:
%  eye = modelEyeParameters()
%
% Description:
%   This routine returns the parameters that define a model eye that is
%   used in the sceneGeometry routines. We make use of values derived by
%   Atchison for human eyes:
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
%   dimensions [depth, horizontal, vertical]; negative values of depth are
%   towards the center of the eye. The model assumes the optical and pupil
%   axis of the eye are algined.
%
% Inputs:
%   none
%
% Optional key/value pairs:
%  'sphericalAmetropia'   - Scalar, in units of diopters. The
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
%  'spectralDomain'       - String, options include {'vis','nir'}.
%                           This is the light domain within which imaging
%                           is being performed. The refractive indices vary
%                           based upon this choice.
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
    eye = modelEyeParameters('sphericalAmetropia',-3,'eyeLaterality','left');
%}


%% input parser
p = inputParser; p.KeepUnmatched = true;

% Optional
p.addParameter('sphericalAmetropia',0,@isnumeric);
p.addParameter('axialLength',[],@(x)(isempty(x) || isnumeric(x)));
p.addParameter('kappaAngle',[],@(x)(isempty(x) || isnumeric(x)));
p.addParameter('eyeLaterality','Right',@ischar);
p.addParameter('species','Human',@ischar);
p.addParameter('spectralDomain','nir',@ischar);

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
        
        
        %% Cornea front surface
        % The properties of the cornea are typically described by the
        % radius of curvature (R) at the vertex and its asphericity (Q).
        % These values are taken from Atchison 2006, Table 1. The radius of
        % curvature of the front surface at the apex varies by spherical
        % ametropia of the eye; Q does not vary.
        eye.corneaFrontSurfaceR = 7.77 + 0.022 * p.Results.sphericalAmetropia;
        eye.corneaFrontSurfaceQ = -0.15;
        
        % We model the cornea as a prolate ellipsoid that is radially
        % symmetric about the optical axis, and require the radii of the
        % ellipsoid. The major and minor radii of an ellipse (a,
        % b) are related to R and Q by:
        %   R = a^2/b
        %	Q = (a^2 / b^2) - 1
        % Therefore, given R and Q, we can obtain a and b, which correspond
        % to the radii of the ellipsoid model, with a corresponding to the
        % axial dimension, and b to the horizontal and verical dimensions.
        % Checking my algebra here:
        %{
            syms a b R Q
            eqn1 = R == a^2/b;
            eqn2 = Q == (a^2 / b^2) - 1;
            solution = solve([eqn1, eqn2]);
            solution.a
            solution.b
        %}
        a = eye.corneaFrontSurfaceR * sqrt( 1 / (eye.corneaFrontSurfaceQ + 1 ) );
        b = eye.corneaFrontSurfaceR / (eye.corneaFrontSurfaceQ + 1 );
        eye.corneaFrontSurfaceRadii(1) = a;
        eye.corneaFrontSurfaceRadii(2:3) = b;
        
        % We set the axial apex of the corneal front surface at position
        % [0, 0, 0]
        eye.corneaFrontSurfaceCenter = [-eye.corneaFrontSurfaceRadii(1) 0 0];
        
        
        %% Cornea back surface
        % The radius of curvature for the back corneal surface was not
        % found to vary by spherical ametropia. The asphericity Q for the
        % back corneal surface was set by Atchison to -0.275.
        eye.corneaBackSurfaceR = 6.4;
        eye.corneaBackSurfaceQ = -0.275;
        
        % Compute the radii of the ellipsoid
        a = eye.corneaBackSurfaceR * sqrt( 1 / (eye.corneaBackSurfaceQ + 1 ) );
        b = eye.corneaBackSurfaceR / (eye.corneaBackSurfaceQ + 1 );
        eye.corneaBackSurfaceRadii(1) = a;
        eye.corneaBackSurfaceRadii(2:3) = b;
        
        % The center of the cornea circle for the back surface is
        % positioned so that there is 0.55 mm of corneal thickness between
        % the front and back surface of the cornea at the apex, following
        % Atchison 2006.
        eye.corneaBackSurfaceCenter = [-eye.corneaFrontSurfaceRadii(1)+0.55 0 0];
        
        
        %% Pupil
        % We position the pupil plane at the depth of the anterior point of
        % the lens. The coordinate space of the model eye is define with
        % respect to the center of the pupil, so the p2 and p3 values are
        % zero
        eye.pupilCenter = [-3.7 0 0];
        
        
        %% Iris
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
        eye.irisRadius = 5.92;
        
        % The iris center is shifted slightly temporally and upward with
        % respect to the pupil center:
        %
        %   ...the typical entrance pupil is decentered
        %   approximately 0.15 mm nasally and 0.1 mm inferior to the
        %   geometric center of the visible iris circumference
        %
        % Bennett, Edward S., and Barry A. Weissman, eds. Clinical contact
        % lens practice. Lippincott Williams & Wilkins, 2005, p119
        %
        % We model an eye with zero iris angle, and thus set the depth of
        % the iris plane equal to the pupil plane.
        switch eyeLaterality
            case 'Right'
                eye.irisCenter = [-3.7 -0.15 -0.1];
            case 'Left'
                eye.irisCenter = [-3.7 0.15 -0.1];
        end
        
        
        %% Posterior chamber
        % The posterior chamber of the eye is modeled as an ellipsoid.
        % Atchison 2005 provides the radii of an ellipsoid model for the
        % posterior chamber and how these dimensions vary with spherical
        % refractive error.
        eye.posteriorChamberRadii = [...
            10.148 - 0.163 * p.Results.sphericalAmetropia ...
            11.365 - 0.090 * p.Results.sphericalAmetropia ...
            11.455 - 0.043 * p.Results.sphericalAmetropia ];
        
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
            %
            % GKA to follow up: Axial length is usually measured with the
            % IOL master along the visual (as opposed to optic or
            % pupillary) axis of the eye. May want to correct for this
            % somewhere.
            scaleFactor = (p.Results.axialLength - posteriorChamberApexDepth) / (eye.posteriorChamberRadii(1)*2);
            eye.posteriorChamberRadii = eye.posteriorChamberRadii .* scaleFactor;
            eye.axialLength = p.Results.axialLength;
        end
        
        % Set the depth of the center of the posterior chamber
        eye.posteriorChamberCenter = ...
            [(-posteriorChamberApexDepth - eye.posteriorChamberRadii(1)) 0 0];
        

        %% Rotation center
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
        
        
        %% Kappa
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
        %	Mathur, Ankit, Julia Gehrmann, and David A. Atchison. "Pupil
        %	shape as viewed along the horizontal visual field." Journal of
        %	vision 13.6 (2013): 3-3.
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
        % Tabernero 2007 report a mean horizontal kappa of 5 degrees in
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
        
        
        %% Refractive indices
        % Obtain refractive index values for this spectral domain.
        eye.corneaRefractiveIndex = returnRefractiveIndex( 'cornea', p.Results.spectralDomain );
        eye.aqueousRefractiveIndex = returnRefractiveIndex( 'aqueous', p.Results.spectralDomain );
        
    otherwise
        error('Please specify a valid species for the eye model');
end

% Meta data regarding the units of the model
eye.meta = p.Results;
eye.meta.units = 'mm';
eye.meta.coordinates = 'eyeWorld';
eye.meta.dimensions = {'depth (axial)' 'horizontal' 'vertical'};
eye.meta.kappa = 'Degrees angle of visual axis w.r.t. pupil axis.';

end % function


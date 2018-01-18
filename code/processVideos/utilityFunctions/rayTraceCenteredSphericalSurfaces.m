function [outputRay, thetas, imageCoords] = rayTraceCenteredSphericalSurfaces(coordsInitial, thetaInitial, opticalSystem, figureFlag)
% Returns the position and angle of a resultant ray WRT optical axis
%
% Syntax:
%   [thetaOut, imagePositionCurrent] = rayTraceCenteredSphericalSurfaces(thetaInitial, nInitial, opticalSystem)

% Description:
%   This routine implements the 2D generalized ray tracing equations of:
%
%       Elagha, Hassan A. "Generalized formulas for ray-tracing and
%       longitudinal spherical aberration." JOSA A 34.3 (2017): 335-343.
%
%   The equations assume a set of spherical surfaces, with each spherical
%   surface having its center of curvature positioned on the optical axis.
%   A ray is specified as originating from position 0 on the optical axis,
%   traveling left-to-right, and making the angle theta with the optical
%   axis. By convention, the optical axis is termed "z", and the orthogonal
%   axis is termed "height". Positive values of theta correspond to the ray
%   diverging to a position above the optical axis. Each spherical surface
%   is specified by a center of curvature and a radius. The center of
%   curvature must lie on the optical axis; positive values place the
%   center to the right of the origin of the ray. A positive radius
%   presents the ray with a convex surface; a negative radius presents the
%   ray with a concave surface. The output of the routine is the position
%   and angle at which the ray (or its reverse projection) intersects the
%   optical axis.
%
% Inputs:
%   coordsInitial         - a 2x1 matrix, with the values corresponding to
%                           the z-position and height of the initial
%                           position of the ray.
%   thetaInitial          - A scalar in radians. A value of zero is aligned
%                           with the optical axis. Values between 0 and pi
%                           direct the ray to diverge "upwards" away from
%                           the axis.
%   opticalSystem         - An mx3 matrix, where m is the number of
%                           surfaces in the model, including the initial
%                           position of the ray. Each row contains the
%                           values [center, radius, refractiveIndex] that
%                           define a spherical lens. The first row
%                           corresponds to the initial conditions of the
%                           ray. Thus, the refractive index value given in
%                           the first row species the index of the medium
%                           in which the raw arises. The center and radius
%                           values for the first row are ignored.
%   figureFlag            - Logical or structure. If logical, the true or
%                           false value controls whether the default style
%                           plot is created. If empty, it will be set to
%                           false. More control of the plot can be obtained
%                           by passing a structure with individual elements
%                           set to true or false.
%
% Outputs:
%   outputRay
%   thetaOut              - A scalar in radians
%   imagePosition         - The point at which the resultant ray (or its
%                           virtual extension) intersects the optical axis.
%
% Examples:
%   Examples 2 and 3 require the modelEyeParametersFunction from the
%   transparentTrack toolbox.
%
%   Ex. 1 - Elagha 2017
%       The paper provides a numerical example in section C which is
%       implemented here as an example. Compare the returned theta values
%       with those given on page 340, section C.
%{
    thetaInitial = deg2rad(17.309724);
    coordsInitial = [0 0];
    figureFlag=true;
    opticalSystem=[nan nan 1; 22 10 1.2; 9 -8 1; 34 12 1.5; 20 -10 1.0];
    [outputRay, thetas, imageCoords] = rayTraceCenteredSphericalSurfaces(coordsInitial, thetaInitial, opticalSystem, figureFlag);
    for ii=1:length(thetas)
        fprintf('theta%d: %f \n',ii-1,rad2deg(thetas(ii)));
    end
    fprintf('Elegha gives a final image distance of 17.768432. we obtain:\n')
    fprintf('i5 - c5 = K4 = %f \n',imageCoords(end,1)-opticalSystem(5,1));
%}
%
%   Ex.2 - Pupil through cornea -
%       A model of the passage of a point on the pupil perimeter through
%       the cornea (units in mm).
%{
    eye = modelEyeParameters();
    pupilRadius = 2;
    theta = deg2rad(-45);
    coords = [eye.pupilCenter(1) pupilRadius];
    opticalSystem = [nan nan eye.aqueousRefractiveIndex; ...
                     eye.corneaBackSurfaceCenter(1) -eye.corneaBackSurfaceRadius eye.corneaRefractiveIndex; ...
                     eye.corneaFrontSurfaceCenter(1) -eye.corneaFrontSurfaceRadius 1.0];
    figureFlag=true;
    rayTraceCenteredSphericalSurfaces(coords, theta, opticalSystem, figureFlag);
%}
%
%   Ex.3 - Pupil through cornea, multiple points and rays -
%{
    eye = modelEyeParameters();
    pupilRadius = 2;
    opticalSystem = [nan nan eye.aqueousRefractiveIndex; ...
                     eye.corneaBackSurfaceCenter(1) -eye.corneaBackSurfaceRadius eye.corneaRefractiveIndex; ...
                     eye.corneaFrontSurfaceCenter(1) -eye.corneaFrontSurfaceRadius 1.0];
    figure;
    clear figureFlag
    figureFlag.show = true;
    figureFlag.new = false;
    figureFlag.surfaces = true;
    figureFlag.imageLines = true;
    figureFlag.rayLines = true;
    figureFlag.textLabels = false;
    for theta = -35:70:35
        for pupilRadius = -2:4:2
            rayTraceCenteredSphericalSurfaces([eye.pupilCenter(1) pupilRadius], theta, opticalSystem, figureFlag);
        end
    end
%}


%% Check the input
% If a value was not passed for figureFlag, set to false.
if nargin==3
    figureFlag.show = false;
    figureFlag.new = false;
    figureFlag.surfaces = false;
    figureFlag.imageLines = false;
    figureFlag.rayLines = false;
    figureFlag.textLabels = false;
end

if nargin==4
    if islogical(figureFlag)
        if figureFlag
            clear figureFlag
            figureFlag.show = true;
            figureFlag.new = true;
            figureFlag.surfaces = true;
            figureFlag.imageLines = true;
            figureFlag.rayLines = true;
            figureFlag.textLabels = true;
        else
            clear figureFlag
            figureFlag.show = false;
            figureFlag.new = false;
            figureFlag.surfaces = false;
            figureFlag.imageLines = false;
            figureFlag.rayLines = false;
            figureFlag.textLabels = false;
        end
    end
end

%% Initialize variables and plotting
nSurfaces = size(opticalSystem,1);
intersectionCoords=zeros(nSurfaces+1,2);
imageCoords=zeros(nSurfaces,2);
thetas=zeros(nSurfaces,1);
relativeIndices=ones(nSurfaces,1);
aVals=ones(nSurfaces,1);

% Set the values for at the first surface (initial position of ray)
thetas(1)=thetaInitial;
intersectionCoords(1,:)=coordsInitial;
imageCoords(1,:)=[coordsInitial(1)-(coordsInitial(2)/tan(thetaInitial)) 0];
opticalSystem(1,1)=imageCoords(1,1);
opticalSystem(1,2)=0;

% Start the figure
if figureFlag.show
    if figureFlag.new
        figure
    end
    hold on
    refline(0,0)
    axis equal
end


%% Peform the ray trace
for ii = 2:nSurfaces
    % The distance between the center of curvature of the current lens
    % surface and the center of curvature of the prior lens surface
    d = opticalSystem(ii,1)-opticalSystem(ii-1,1);
    % The relative refractive index of the prior medium to the medium of
    % the surface that the ray is now impacting
    relativeIndices(ii)=opticalSystem(ii-1,3)/opticalSystem(ii,3);
    % implements equation 54 of Eleghan
    aVals(ii) = ...
        (1/opticalSystem(ii,2))*(relativeIndices(ii-1).*aVals(ii-1).*opticalSystem(ii-1,2)+d.*sin(thetas(ii-1)));
    % check if the incidence angle is above the critical angle for the
    % relative refractive index at the surface interface
    if abs((aVals(ii)*relativeIndices(ii))) > 1
        warning('Angle of incidence for surface %d greater than critical angle',ii);
        break
    end
    % Find the angle of the ray after it enters the current surface
    thetas(ii) = thetas(ii-1) - asin(aVals(ii)) + asin(aVals(ii).*relativeIndices(ii));
    % Find the coordinates at which the ray, after making contact with the
    % current surface, would contact (or originate from) the optical axis
    imageCoords(ii,:) = [opticalSystem(ii,1) + ...
        1./(-(1/relativeIndices(ii)).*(1/(aVals(ii).*opticalSystem(ii,2))).*sin(thetas(ii))) 0];
    % Find the coordinate of at which the ray intersected the current surface
    if ii==2
        slope = tan(thetas(ii-1));
    else
        slope = intersectionCoords(ii-1,2)/(intersectionCoords(ii-1,1)-imageCoords(ii-1,1));
    end
    intercept = (0-imageCoords(ii-1))*slope;
    % the linecirc routine returns coordinates of intersection
    % between the ray and a circle. I
    [xout,yout] = linecirc(slope,intercept,opticalSystem(ii,1),0,abs(opticalSystem(ii,2)));
    % This next bit of logic figures out which of the two coordinates of
    % intersection of the ray with a sphere correspond to the one we want
    % for the lens
    if length(xout)==2
        whichIdx = 1.5+0.5*((sign(opticalSystem(ii,1))*sign(opticalSystem(ii,2))));
        intersectionCoords(ii,:)=[xout(whichIdx) yout(whichIdx)];
    else
        % If it returns only one coordinate the ray either was tangential to
        % the surface or missed entirely. We therefore exit the ray tracing
        break
    end
    % Update the plot
    if figureFlag.show
        % add this lens surface
        if figureFlag.surfaces
            plotLensArc(opticalSystem(ii,:))
        end
        % plot the line for the virtual image
        if figureFlag.imageLines
            plot([imageCoords(ii-1,1) intersectionCoords(ii-1,1)],[imageCoords(ii-1,2) intersectionCoords(ii-1,2)],'--b');
        end
        % plot the line for the path of the ray
        if figureFlag.rayLines
            plot([intersectionCoords(ii-1,1) intersectionCoords(ii,1)],[intersectionCoords(ii-1,2) intersectionCoords(ii,2)],'-r');
        end
    end
end
% Extend the final ray as a unit vector
slope = tan(thetas(ii));
intersectionCoords(ii+1,:)=[intersectionCoords(ii,1)+1 intersectionCoords(ii,2)+slope];


%% Finish and clean up
% Assemble an output which is the unit vector for the final ray
outputRay = [intersectionCoords(ii,:); intersectionCoords(ii+1,:)];

% Complete the plot
if figureFlag.show
    % Plot the final virtual image
    if figureFlag.imageLines
        plot([imageCoords(ii,1) intersectionCoords(ii,1)],[imageCoords(ii,2) intersectionCoords(ii,2)],'-b');
    end
    % Plot the output unit ray vector
    if figureFlag.rayLines
        plot([intersectionCoords(ii,1) intersectionCoords(ii+1,1)],[intersectionCoords(ii,2) intersectionCoords(ii+1,2)],'-r');
    end
    % Replot the refline    
    refline(0,0)
    % Add some labels
    if figureFlag.textLabels
        plot(opticalSystem(:,1),zeros(nSurfaces,1),'+k');
        text(opticalSystem(:,1),zeros(nSurfaces,1)-(diff(ylim)/50),strseq('c',[1:1:nSurfaces]),'HorizontalAlignment','center')
        plot(imageCoords(:,1),zeros(nSurfaces,1),'*r');
        text(imageCoords(:,1),zeros(nSurfaces,1)-(diff(ylim)/50),strseq('i',[1:1:nSurfaces]),'HorizontalAlignment','center')
    end
    hold off
end

end % function


%% LOCAL FUNCTIONS
function plotLensArc(opticalSystem)
ang=pi/2:0.01:3*pi/2;
xp=opticalSystem(2)*cos(ang);
yp=opticalSystem(2)*sin(ang);
plot(opticalSystem(1)+xp,yp,'-k');
end
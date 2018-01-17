function [thetaOut, imagePositionCurrent] = rayTraceCenteredSphericalSurfaces(thetaInitial, nInitial, opticalSystem)
% Returns the position and angle of a resultant ray WRT optical axis
%
% Syntax:
%   [thetaOut, imagePositionCurrent] = rayTraceCenteredSphericalSurfaces(thetaInitial, nInitial, opticalSystem)

% Description:
%   This routine implements the generalized ray tracing equation of:
%
%       Elagha, Hassan A. "Generalized formulas for ray-tracing and
%       longitudinal spherical aberration." JOSA A 34.3 (2017): 335-343.
%
%   The equations assume a set of spherical surfaces, with each spherical
%   surface having its center of curvature positioned on the optical axis.
%   A ray is specified as originating from position 0 on the optical axis,
%   traveling left-to-right, and making the angle theta with the optical
%   axis. Positive values of theta correspond to the ray diverging to a
%   position above the optical axis. Each spherical surface is specified by
%   a center of curvature and a radius. The center of curvature must lie on
%   the optical axis; positive values place the center to the right of the
%   origin of the ray. A positive radius presents the ray with a convex
%   surface; a negative radius presents the ray with a concave surface.
%   The output of the routine is the position and angle at which the ray
%   (or its reverse projection) intersects the optical axis.
%
% Inputs:
%   thetaInitial          - A scalar in radians. A value of zero is aligned
%                           with the optical axis. Values between 0 and pi
%                           direct the ray to diverge "upwards" away from
%                           the axis.
%   nInitial              - The refractive index of the medium in which the
%                           ray arises.
%   opticalSystem         - An mx3 matrix. Each row contains the values
%                           [center, radius, refractiveIndex] that define a spherical lens.
%
% Outputs:
%   thetaOut              - A scalar in radians
%   imagePosition         - The point at which the resultant ray (or its
%                           virtual extension) intersects the optical axis.
%
% Examples:
%   Elagha provides a numerical example in section C of the paper. This
%   validation code replicates the example. The returned values should be
%   thetaOut = -26.583586 degrees, and imagePosition = 17.768432.
%{
    thetaInitial = deg2rad(17.309724);
    nInitial=1;
    opticalSystem=[22,10,1.2; 9,-8,1; 34,12,1.5; 20,-10,1];
    [ thetaOut, imagePosition ] = rayTraceCenteredSphericalSurfaces(thetaInitial, nInitial, opticalSystem );
    rad2deg(thetaOut), imagePosition
%}

figure
refline(0,0)
hold on

aPrior = [];
imagePositionPrior = [0 0];
intersectPointPrior = [0 0];
thetaIn = thetaInitial;
thetaOut = [];

nLenses = size(opticalSystem,1);
for ii = 1:nLenses
    if ii==1
        aCurrent = (opticalSystem(ii,1).*sin(thetaIn))/opticalSystem(ii,2);
        nRelCurrent = nInitial/opticalSystem(ii,3);
    else
        % distance between the center of the current lens and the center of
        % the prior lens
        d = opticalSystem(ii,1)-opticalSystem(ii-1,1);
        % relative refractive index of the prior medium to the medium in
        % which the incident ray is currently located
        if ii==2
            nRelPrior = nInitial/opticalSystem(ii-1,3);
        else
            nRelPrior = opticalSystem(ii-2,3)/opticalSystem(ii-1,3);
        end
        % relative refractive index of the medium in which the incident ray
        % is currently located to the medium of the lens it is about to
        % intersect
        nRelCurrent = opticalSystem(ii-1,3)/opticalSystem(ii,3);
        
        aCurrent = (1/opticalSystem(ii,2))*(nRelPrior.*aPrior.*opticalSystem(ii-1,2)+d.*sin(thetaIn));
    end
    % Find the angle of the ray after it enters the current lens
    thetaOut = thetaIn - asin(aCurrent) + asin(aCurrent.*nRelCurrent);
    % Find the position at which the ray (or its virtual extension) would
    % contact the optical axis
    imagePositionCurrent = [1./(-(1/nRelCurrent).*(1/(aCurrent.*opticalSystem(ii,2))).*sin(thetaOut)) 0];
    % Find the position of the ray at its point of contact on the current
    % lens
    z = ((imagePositionCurrent(1)-imagePositionPrior(1)).*sin(thetaIn))./sin(pi-thetaIn+thetaOut);
    intersectPointCurrent=[imagePositionCurrent-cos(-thetaOut).*z sin(-thetaOut).*z];
    % Update the plot
    plot([intersectPointPrior(1) intersectPointCurrent(1)],[intersectPointPrior(2) intersectPointCurrent(2)],'-k');
    if ii==nLenses
    else
        thetaIn=thetaOut;
        aPrior=aCurrent;
        imagePositionPrior = imagePositionCurrent;
        intersectPointPrior = intersectPointCurrent;
    end
end
plot([intersectPointCurrent(1) imagePositionCurrent(1)],[intersectPointCurrent(2) imagePositionCurrent(2)],'-r');
end



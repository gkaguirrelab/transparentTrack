function [ecc, pol] = convertScreenCoordinatesToPolar(x,y,f)
% [ecc, pol] = convertScreenToPolar(x,y,f)
%
% this function converts 3-D screen coordinates to 3-D polar coordinates.
% The coordinate systems share the origin (center of the screen) and the
% Z coordinates (viewing distance, measured from the center of the screen).
% 
% IMPORTANT NOTE: matlab's "polar" functions assume a different orientation
% for the cartesian and the polar plots. Pay attention when plotting these
% results:
% - flip the Y axis upsidedown when plotting screen coordinates 
% - flip the polar plot to be clockwise and orient the zero on the upper
% vertical meridian
%
% Inputs
% x = X coordinate of the gaze (orizontal growing from left to right)
% y = Y coordinate of the gaze (vertical growing from top to bottom)
% f = viewing distance
% Note that all inputs need to be in the same units (e.g. mm)
%
% Outputs
% ecc = Polar radius in degrees of visual angle of the  gaze
% direction from the center of the screen.
% pol = polar position of the gaze, given as clockwise from 0
% degrees (with the upper vertical meridian being the 0 degree position,
% and the right, horizontal meridian being 90 degrees).
% 


if isnan(x) || isnan(y)
    ecc = NaN;
    pol = NaN;

else
    xa = abs(x);
    d = hypot(x,y);
    
    ecc = atand(d/(f));
    
    if x==0 && y==0
        pol = 0;
        
    elseif x==0 && y>0
        pol = 180;
        
    elseif x==0 && y<0
        pol = 0;
        
    elseif y==0 && x>0
        pol = 90;
        
    elseif y==0 && x<0
        pol = 270;
        
    elseif x>0 && y>0
        deltaAng = 180;
        pol = deltaAng - asind(xa/d);
        
    elseif x>0 && y<0
        deltaAng = 0;
        pol = deltaAng + asind(xa/d);
        
    elseif x<0 && y>0
        deltaAng = 180;
        pol = deltaAng + asind(xa/d);
        
    elseif x<0 && y<0
        deltaAng = 360;
        pol = deltaAng - asind(xa/d);
        
    end
end

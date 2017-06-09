function [ecc, pol] = liveTrackCartToPol(x,y,f)
% function [ecc, pol] = ConvertLiveTrackCartToPolar(x,y,f)
%
% this function converts the cartesian coordinates of the Gaze direction
% from the LiveTrack Report reference system to the standard eccentricity
% and polar position values.
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
% 
% 
%  Jan 2017 - written by Giulia Frazzetta 

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

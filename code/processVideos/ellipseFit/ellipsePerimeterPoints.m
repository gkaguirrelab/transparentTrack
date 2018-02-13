function [ Xp, Yp ] = ellipsePerimeterPoints( transparentEllipseParams, steps )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if nargin == 1
    steps = 5;
end

p = ellipse_transparent2ex(transparentEllipseParams);
x = p(1);
y = p(2);
a = p(3);
b = p(4);
theta = p(5);

sintheta = sin(theta);
costheta = cos(theta);

alpha = linspace(0, 2*pi-(2*pi/steps), steps)';
sinalpha = sin(alpha);
cosalpha = cos(alpha);

Xp = x + (a * cosalpha * costheta - b * sinalpha * sintheta);
Yp = y + (a * cosalpha * sintheta + b * sinalpha * costheta);

end


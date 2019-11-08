function [ curve ] = circle(r,c)
% creates the object of class Curve corresponding to the circle C(r,c),
% where r is the radius and c the center

x = @(t)(c(1) + r*cos(t));
y = @(t)(c(2) + r*sin(t));
I = [-pi,pi];
boundedSide = 'left';
curve = SimpleCurve(x,y,I,boundedSide);


end


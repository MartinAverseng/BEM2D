function [ curve ] = openline(a,b)

x = @(t)(t);
y = @(t)(0*t);
I = [a,b];
curve = SimpleCurve(x,y,I);

end


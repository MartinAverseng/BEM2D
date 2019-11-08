function[c] = boomerang()

x = @(t)((0.7 + cos(t*pi)).*cos(t*pi));
y = @(t)(sin(t*pi));
c = SimpleCurve(x,y,[-1,1],'left');

end
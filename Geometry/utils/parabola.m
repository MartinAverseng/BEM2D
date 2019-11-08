function[c,incWave] = parabola(k)

x = @(t)(t);
y = @(t)(t.^2 - 0.25);
I = [-1,1];
c = SimpleCurve(x,y,I);

if and(nargout == 2,nargin ==1)
    theta_inc = pi/2;
    X = R2toRfunc.X; Y = R2toRfunc.Y;
    incWave = exp(1i*k*(X*cos(theta_inc) + Y*sin(theta_inc)));
end

end
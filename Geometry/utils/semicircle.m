function[c,incWave,dxf,dyf] = semicircle(k)

x = @(t)(cos(pi/2*t));
y = @(t)(sin(pi/2*t));
I = [-1, 1];
c = SimpleCurve(x,y,I);

if and(nargout == 2,nargin ==1)
    theta_inc = 0;
    X = R2toRfunc.X; Y = R2toRfunc.Y;
    incWave = exp(1i*k*(X*cos(theta_inc) + Y*sin(theta_inc)));
end
if and(nargout == 4,nargin ==1)
    theta_inc = 0;
    X = R2toRfunc.X; Y = R2toRfunc.Y;
    incWave = exp(1i*k*(X*cos(theta_inc) + Y*sin(theta_inc)));
    dxf = 1i*k*cos(theta_inc)*incWave;% = d(incwave)/dx
    dyf = 1i*k*sin(theta_inc)*incWave;
end


end
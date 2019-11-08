function [arc,incWave ,dxf,dyf] = archSpirale(k)

m = 3;
y = @(s)((1 + 0.5*m*pi*(s + 1)).*cos(m*pi*(s+1))/5);
x = @(s)(-(1 + 0.5*m*pi*(s + 1)).*sin(m*pi*(s+1))/5);

I = [-1,1];
arc = SimpleCurve(x,y,I);

if and(nargout >= 2,nargin ==1)
    theta_inc = -pi/2 + pi/8;
    X = R2toRfunc.X; Y = R2toRfunc.Y;
    incWave = exp(1i*k*(X*cos(theta_inc) + Y*sin(theta_inc)));
    dxf = 1i*k*cos(theta_inc)*incWave;% = d(incwave)/dx
    dyf = 1i*k*sin(theta_inc)*incWave;
end

end


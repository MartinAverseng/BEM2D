function [arc,incWave ,dxf,dyf] = spirale(k)

a= 1;
b= 0.2;
c= 2;
s0 = -0.2;
x = @(s)(a*exp(b*c*(s+s0)).*cos(c*(s+s0)));
y = @(s)(a*exp(b*c*(s+s0)).*sin(c*(s+s0)));
I = [-1,1];
arc = SimpleCurve(x,y,I);

if and(nargout >= 2,nargin ==1)
    theta_inc = pi/4;
    X = R2toRfunc.X; Y = R2toRfunc.Y;
    incWave = exp(1i*k*(X*cos(theta_inc) + Y*sin(theta_inc)));
    dxf = 1i*k*cos(theta_inc)*incWave;% = d(incwave)/dx
    dyf = 1i*k*sin(theta_inc)*incWave;
end

if nargin == 0
    close all
    plot(arc)
end


end




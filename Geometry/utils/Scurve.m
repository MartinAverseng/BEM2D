function [arc,incWave ] = Scurve(k)

a= 0.2;
b= 0.2;
c= 4;
s0 = 0;
x = @(s)(s);
y = @(s)(a*exp(b*c*(s+s0)).*sin(c*(s+s0)));
I = [-1,1];
arc = SimpleCurve(x,y,I);

if and(nargout == 2,nargin ==1)
    theta_inc = pi/8;
    X = R2toRfunc.X; Y = R2toRfunc.Y;
    incWave = exp(1i*k*(X*cos(theta_inc) + Y*sin(theta_inc)));
end



end


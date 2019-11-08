function [arc,incWave ] = Vcurve(k)


x = @(s)(abs(3*s))/3;
y = @(s)(s)/3;
I = [-1,1];
arc = SimpleCurve(x,y,I);

if and(nargout == 2,nargin ==1)
    theta_inc = pi/12;
    X = R2toRfunc.X; Y = R2toRfunc.Y;
    incWave = exp(1i*k*(X*cos(theta_inc) + Y*sin(theta_inc)));
end



end


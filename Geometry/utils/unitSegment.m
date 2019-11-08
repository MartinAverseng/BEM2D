function [ c , incWave, dxf, dyf] = unitSegment(k)

c = openline(-1,1);
X = R2toRfunc.X;
incWave = 1./sqrt(((X^2 + 0.001^2)));


if and(nargout >= 2,nargin ==1)
    theta_inc = 0;
    X = R2toRfunc.X; Y = R2toRfunc.Y;
    incWave = exp(1i*k*(X*cos(theta_inc) + Y*sin(theta_inc)));%*1/sqrt((X^2 + (Y-0.01)^2));
    
    dxf = 1i*k*cos(theta_inc)*incWave;% = d(incwave)/dx
    dyf = 1i*k*sin(theta_inc)*incWave;
end


end

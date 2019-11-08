function [dn_u] = explicitDtN_planeWave_0inc(k,M)
% 
X = R2toRfunc(@(Z)(Z(:,1)));
Y = R2toRfunc(@(Z)(Z(:,2)));
r = sqrt(X^2 + Y^2);
phi = atan2(Y,X);

dn_u = R2toRfunc;
epsilon = [1;2*ones(M,1)];
for m = 0:M
    em = R2toRfunc(@(Z)(cos(m*atan2(Z(:,2),Z(:,1)))));
    dn_u = dn_u + k*epsilon(m+1)*(-1i)^m*(besselj(m,k)/besselh(m,k)*(besselh(m-1,k)-besselh(m+1,k))/2 - 1/2*(besselj(m-1,k) - besselj(m+1,k)))*em;
end

end


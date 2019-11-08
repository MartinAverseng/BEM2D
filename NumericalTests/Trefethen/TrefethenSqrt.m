function [ v ] = TrefethenSqrt(A,N,u,M,min_m,max_M)
% Approximates v = sqrt(A)*u using trefethen's method.
% N is the number of iterations used.
if ~exist('M','var')||isempty(max_M)
    M = 1;
end

k2 = min_m/max_M; % elliptic functions parameter k^2
Kp = ellipke(1-k2);
t = 1i*(.5:N)*Kp/N;
[sn, cn, dn] = ellipj(imag(t),1-k2);
cn = 1./cn; dn = dn.*cn; sn = 1i*sn.*cn;
w = sqrt(min_m)*sn;
dzdt = cn.*dn;

if ~isempty(u)
    
    v = zeros(size(A,1),1);
    for j = 1:N
        v = v - dzdt(j)*((A-w(j)^2*M)\(M*u));
    end
    v = (-2*Kp*sqrt(min_m)/(pi*N))*(M\(A*v));
    v = M*v;
else
    
    Mat = zeros(size(A));
    for j = 1:N
        Mat = Mat - dzdt(j)*inv(A-w(j)^2*M)*M;
    end
    Mat = (-2*Kp*sqrt(min_m)/(pi*N))*inv(M)*(A*Mat);
    Mat = M*Mat;
    v = Mat;
end


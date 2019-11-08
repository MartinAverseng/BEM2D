clear all;
close all;
tic;
X = [0.5 1];
A = [2 0];
B = [2.0001 0];
lnK = logSingK;
sA = 10;
sB = sA + norm(B-A);
u = (B-A)/norm(B-A);
y1 = @(s)(A(1) + (s-sA)*u(1));
y2 = @(s)(A(2) + (s-sA)*u(2));
r = @(s)(sqrt((X(1)-y1(s)).^2 + (X(2)-y2(s)).^2));
fun = @(s)(log(r(s)));
[xhat1,what1] = gaussQuad(fun,4,sA,sB);
s1 = xhat1;
sref1 = (xhat1 - sA)/(sB-sA);
[xhat2,what2] = lnK.singularKernelQuadrature(X,A,B,4);
s2 = xhat2*(sB-sA) + sA;
sref2 = xhat2; % belongs to [0,1]
toc
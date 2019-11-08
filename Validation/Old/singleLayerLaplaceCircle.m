% Test the new version of Op and BIOGalerkineMain;
clear all
close all

%% Validation (We have explicit formulas for spherical harmonics)
r = 1; % do not change r
theCurve = circle(r,[0,0]);
tol = 1e-3;
kernel = LogKernel(1/2); 

% 1°) Formulas for incident and scattered waves 

X = R2toRfunc(@(Z)Z(:,1));
Y = R2toRfunc(@(Z)Z(:,2));


% Development in spherical harmonics
ms = [2 3 10];
coeffs = [1 -1 5];
u0_func{1} = X^2 - Y^2; % Spherical harmonics 
u0_func{2} = X^3 - 3*X*Y^2;
u0_func{3} = X^10 -Y^10 + 45*(-X^8*Y^2 + Y^8*X^2) +210*(X^6*Y^4 - Y^6*X^4);
lambda_theo = R2toRfunc;
u0 = R2toRfunc;
for i = 1:length(ms)
    m = ms(i);
    u0 = u0 + coeffs(i)*u0_func{i};
    % S phi = 1/(2m) phi when phi is a spherical harmonic of degree m
    lambda_theo = lambda_theo + 2*m*coeffs(i)*u0_func{i};
end

Ns = [15 150 500 250 500 1000];
q = 3;

for N_i = 1:length(Ns)
    N = Ns(N_i);
    mesh = MeshCurve(theCurve,N);
    Xh = FEspace(mesh,'P1',q);
    Vh = FEspace(mesh,'P1',q);
    lambda_theoPih = Xh.Pi_h(lambda_theo);
    l = Vh.secondMember(u0);
    [S,Aop] = BIOGalerkine(Xh,Vh,'U',kernel,'V','tol',tol);
    [Sprecond,AopPrecond] = BIOGalerkine(Xh,Vh,'U',kernel,'V','a_factor',3,'noFarField',true);
    S = -1/(2*pi)*(S + regularize(Xh,Vh,'U','ln','V','refineQuad',10));
    Sprecond = -1/(2*pi)*(Sprecond + regularize(Xh,Vh,'U','ln','V','refineQuad',10));
    M = Sprecond.concretePart;
end








%[Sfull,AopFull] = BIOGalerkine(Xh,Vh,'U',kernel,'V','fullMatrix',true);

% Validation of Aop.
% B = AopFull - Aop;
% E1 = [1; zeros(size(B,2)-1,1)];
% C = B*E1;
% plot(abs(C));
% hold on
% plot(0*abs(C) + tol,'--');

% Assembling the single layer operator. 
%Sfull = -1/(2*pi)*(Sfull + regularize(Xh,Vh,'U','ln','V','refineQuad',6));






%% Test for 



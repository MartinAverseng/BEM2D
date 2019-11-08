%% Test de l'idée d'intégrales singulières:
% On cherhe à calculer \int_{seg} ln|x - y|/omega(y) dy 

a = -1; b = -0.999;
x = (a+3*b)/4;

omega = @(x)(sqrt(1-x.^2));

% Intégrale approchée par Matlab 
IMatlab = integral(@(y)(log(abs(x - y))./omega(y)),a,b);


% Intégrale approchée de François : 
l = sqrt((b-a)^2 + (omega(b) - omega(a))^2);
F = @(x)(x.*log(abs(x)) - x);
theta0 = acos(x);
thetaA = acos(a); thetaB = acos(b);
s0 = (theta0 - thetaA)/(thetaB - thetaA)*l;
xprime = a + s0/l*(b-a);
delta = x - xprime; 
IFrancois1 = l/(b-a)*(F(b-x) - F(a - x));
IFrancois2 = l/(b-a)*(F(b-xprime) - F(a - xprime));

disp(abs(IMatlab - IFrancois1));
disp(abs(IMatlab - IFrancois2));
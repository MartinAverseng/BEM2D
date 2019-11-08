%% Test case on the L-shaped domain.
clear all
close all;

LshapedVertices = [0 0; 1 0; 1 1/2; 1/2 1/2; 1/2 1; 0 1];
LshapedDomain = polygonCurve(LshapedVertices,true,'left');
LshapedDomain = ellipse(1,3);

N = 10000;
mesh = MeshCurve(LshapedDomain,N);
plot(mesh);

tol = 1e-5;
q = 3;
Vh = FEspace(mesh,'P1',q);

k = 200;
theta_inc = 2*pi/3;
X =  R2toRfunc(@(Z)(Z(:,1)));
Y =  R2toRfunc(@(Z)(Z(:,2)));
planeWave = exp(1i*k*(X*cos(theta_inc) + Y*sin(theta_inc)));
dx_planeWave = 1i*k*cos(theta_inc)*exp(1i*k*(X*cos(theta_inc) + Y*sin(theta_inc)));
dy_planeWave = 1i*k*sin(theta_inc)*exp(1i*k*(X*cos(theta_inc) + Y*sin(theta_inc)));
g = Vh.Pi_h(planeWave);
figure
plot(g);
drawnow

M = Vh.Mass;
dM = Vh.dMass;
invM = @(u)(M\u);

uN = Vh.Pi_h_normalDerivative(dx_planeWave,dy_planeWave);

X = Vh.gaussPoints;
Sop = singleLayer(Vh,'Op_opt',{'a_factor',20,'tol',1e-3},'k',k);
Sgalerk= Sop.galerkine(Vh,'U');
l = Vh.secondMember(-planeWave) ;

Prec = @(u)(invM(Sgalerk*invM((dM - k^2*M)*invM(u))));


[lambda1,FLAG1,RELRES1,ITER1,RESVEC1] = variationalSol(Sgalerk + 1/2*M,l,[],1e-8,50);
[lambda2,FLAG2,RELRES2,ITER2,RESVEC2] = variationalSol(Sgalerk,l,[],1e-8,50);
figure
semilogy(1:length(RESVEC1),RESVEC1,'-o');
hold on
semilogy(1:length(RESVEC2),RESVEC2,'-o');
drawnow;

xlabel('Nombre d''iterations')
ylabel('Erreur résiduelle')
legend({'Avec Préconditionneur','Sans préconditionneur'})
title('Résolution du système linéaire par méthode itérative')

%% Compute the radiating solution. 

x1 = linspace(-3,3,500)+1/2;
x2 = linspace(-3,3,500)+1/2;

[X1,X2] = meshgrid(x1,x2);
Sop.set_X([X1(:) X2(:)]);
valsDiffr = Sop*(lambda2);
valsInc = planeWave([X1(:) X2(:)]);
valsDiffr = reshape(valsDiffr,size(X1,1),size(X1,2));
valsInc = reshape(valsInc,size(X1,1),size(X1,2));



%% Animate the wave

amplitude = valsInc+valsDiffr;

figure
imagesc(x1,x2,20*log(abs(amplitude)));
axis xy;
axis equal

animateWave(x1,x2,k,amplitude)


%% 


Main;
N = 500;
q = 3;
curve = circle(1,[0,0]);
mesh = MeshCurve(curve,N);
tol = 1e-5;
Vh = FEspace(mesh,'P1',q);

k = 17;

X =  R2toRfunc(@(Z)(Z(:,1)));
Y =  R2toRfunc(@(Z)(Z(:,2)));
planeWave = exp(-1i*k*X);
g = Vh.Pi_h(planeWave);
figure
g.show;
drawnow


X = Vh.gaussPoints;
Sop = singleLayer(k,Vh,Vh.gaussPoints,'a_factor',5,'tol',1e-3);
Sgalerk = Sop.galerkine(Vh,'U');
l = Vh.secondMember(-planeWave);

[lambda,FLAG,RELRES,ITER,RESVEC] = variationalSol(Sgalerk,l,[],1e-8,300);
figure
semilogy(1:length(RESVEC),RESVEC,'-o');
drawnow;

figure
show(real(lambda));
hold on;
dn_u = explicitDtN_planeWave_0inc(k,150);
lambda_EXACT = Vh.Pi_h(dn_u);
show(real(lambda_EXACT))
drawnow;

%% Compute the radiating solution. 

x1 = linspace(-3,3,5000);
x2 = linspace(-3,3,500);

[X1,X2] = meshgrid(x1,x2);
R2 = X1.^2 + X2.^2;
Sop.set_X([X1(:),X2(:)],'a_factor',1/4);
valsDiffr = Sop*lambda;
valsDiffr(R2<1)=0;
valsInc = planeWave([X1(:) X2(:)]);
valsInc(R2<1) = 0;
valsDiffr = reshape(valsDiffr,size(X1,1),size(X1,2));
valsInc = reshape(valsInc,size(X1,1),size(X1,2));

%% Animate the wave

amplitude = valsInc+valsDiffr;
animateWave(x1,x2,k,amplitude)



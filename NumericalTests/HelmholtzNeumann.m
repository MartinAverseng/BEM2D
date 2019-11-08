%% Test case on the some shape with edges.

Main;

nn = 50;
k = nn*pi/2;
[curve,incWave,dxf,dyf] = unitSegment(k);

figure
plot(curve);
l = curve.length;
N = fix(10*k);
meshAdapt = MeshCurve(curve,N,@cos,[-pi,0]);
Vh =  weightedFEspace(meshAdapt,'P1','1/sqrt(1-t^2)',...
    'quadNum',3,'specialQuadSegs',1:meshAdapt.nseg);

Wh =  weightedFEspace(meshAdapt,'P1','sqrt(1-t^2)',3,'specialQuadSegs',1:meshAdapt.nseg);
M = Wh.Mass.concretePart;
[L,U,P,Q] = lu(M);
invM = @(u)(Q*(U\(L \(P*u))));

dM = (Vh.omega_dx_omega)'*AbstractMatrix.spdiag(Vh.W)*(Vh.omega_dx_omega);
dM = dM.concretePart;
omega2 = Wh.omega2;
K1 = dM - k^2*(omega2-M);
K = dM - k^2*omega2;
Np = 15;
theta = pi/3;
keps = k+1i*0.025*k^(1/3);
sqrtDarbasK1 = @(x)(padePrecondDarbas(x,Np,theta,keps,M,K1));
clear dM;
clear omega2;

Op_opt = {'tol',1e-3,'fullMatrix',true};
Nw = hyperSingular_w(Vh,'k',k,'Op_opt',Op_opt);
Nwgalerk = Nw.galerkine;

PrecDarbas = @(u)(K\sqrtDarbasK1(M\u));

clear M L U Q P dM K K1 omega2 
secondMemb = Wh.normalDerivative(dxf,dyf);

t0 = tic;
[lambda0,FLAG0,RELRES0,ITER0,RESVEC0] = variationalSol(Nwgalerk,secondMemb,[],1e-8,N);
t0 = toc(t0);
disp(t0)

t1 = tic;
[lambda1,FLAG1,RELRES1,ITER1,RESVEC1] = variationalSol(Nwgalerk,secondMemb,[],1e-8,N,PrecDarbas);
t1 = toc(t1);
disp(t1);


figure
semilogy(1:length(RESVEC0),RESVEC0/norm(secondMemb.concretePart),'-o');
hold on
semilogy(1:length(RESVEC1),RESVEC1/norm(PrecDarbas(secondMemb.concretePart)),'--x');


xlabel('Iteration number')
ylabel('Residual error')

legend({'Without preconditioner','With preconditioner'});
legend boxoff

set(gca,'FontSize',15)

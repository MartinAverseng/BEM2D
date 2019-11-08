%% Test case on the some shape with edges.

Main;
k = 0;
curve = unitSegment();
X = R2toRfunc.X;
incWave = sqrt(((X^2 + 0.001^2)));

figure
plot(curve);
l = curve.length;
N = 1600;
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

Op_opt = {'fullMatrix',true};

Nw = hyperSingular_w(Vh,'k',k,'Op_opt',Op_opt);
Nwgalerk = Nw.galerkine;

PrecTref = @(u)(dM\TrefethenSqrt(dM,15,invM(u),M,1/2,Vh.ndof^2));
clear M L U Q P dM  
secondMemb = Wh.secondMember(incWave);

t0 = tic;
[lambda0,FLAG0,RELRES0,ITER0,RESVEC0] = variationalSol(Nwgalerk,secondMemb,[],1e-8,N);
t0 = toc(t0);
disp(t0)

t1 = tic;
[lambda1,FLAG1,RELRES1,ITER1,RESVEC1] = variationalSol(Nwgalerk,secondMemb,[],1e-8,N,PrecTref);
t1 = toc(t1);
disp(t1);

figure
semilogy(1:length(RESVEC0),RESVEC0/norm(secondMemb.concretePart),'-o');
hold on
semilogy(1:length(RESVEC1),RESVEC1/norm(PrecTref(secondMemb.concretePart)),'--x');


xlabel('Iteration number')
ylabel('Residual error')

legend({'Without preconditioner','With preconditioner'})
legend boxoff

set(gca,'FontSize',15)
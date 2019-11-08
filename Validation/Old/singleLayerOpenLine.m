% Validate the single layer potential on open line segment.

Main;
N = 200;
n = 0;
q = 30;
curve = openline(-1,1);
mesh = MeshCurve(curve,N);

Xh = FEspace(mesh,'P0',q);
Vh = FEspace(mesh,'P1',q);
Wh = FEspace(mesh,'P2',q);


S0 = Xh.singleLayerPotential();
S1 = Vh.singleLayerPotential();
S2 = Wh.singleLayerPotential();

if n==0
    lambda_n = log(2)/2;
else
    
    lambda_n = 1/(2*n);
end

u0 = R2toRfunc(@(Z)(chebyshevT(n,Z(:,1))));
lambda_theo = R2toRfunc(@(x)(1./(lambda_n*sqrt(1.0000000000000001 - x(:,1)).*sqrt(1.0000000000000001 +x(:,1))).*chebyshevT(n,x(:,1))));

l0 = Xh.secondMember(u0);
l1 = Vh.secondMember(u0);
l2 = Wh.secondMember(u0);

lambda_theo_FE = func(Wh,lambda_theo);
t0 = tic;
lambda0 = variationalSol(S0,l0);
t0 = toc(t0);
t1 = tic;
lambda1 = variationalSol(S1,l1);
t1 = toc(t1);
t2 = tic;
lambda2 = variationalSol(S2,l2);
t2 = toc(t2);


lambda0.show;
hold on
lambda1.show;
lambda2.show;
lambda_theo_FE.show;
xlim([0.1 1.9]);
ylim('auto');






















% Validate the single layer potential on the circle.











%%
Main;
N = 500;
q = 3;
r = 1; % do not change r
curve = circle(r,[0,0]);
mesh = MeshCurve(curve,N);
tol = 1e-6;

Xh = FEspace(mesh,'P0',q);
Vh = FEspace(mesh,'P1',q);
Wh = FEspace(mesh,'P2',q);

k = 30;

S0k = singleLayer(k,Xh,Xh,'refineQuad',15,'tol',tol);
S1k = singleLayer(k,Vh,Vh,'refineQuad',15,'tol',tol);
M0k = S0k.concretePart;
M1k = S1k.concretePart;

m_list = [3 10 13];
m_coeff = [1 1 2];

u0 = R2toRfunc;
lambda_theo = R2toRfunc;
for mi = 1:length(m_list)
    m = m_list(mi);
    em = R2toRfunc(@(Z)(cos(m*atan2(Z(:,2),Z(:,1)))));
    u0 = u0 + m_coeff(mi)*besselj(m,k)*em;
    lambda_theo = lambda_theo + k*m_coeff(mi)*(-besselj(m+1,k)...
        +besselh(m+1,1,k)*besselj(m,k)/besselh(m,1,k))*em;
end



l0 = Xh.secondMember(u0);
l1 = Vh.secondMember(u0);
l2 = Wh.secondMember(u0);

t0 = tic;
lambda0 = S0k\l0;
lambda0theo_Pih = Xh.Pi_h(lambda_theo);

t0 = toc(t0);
t1 = tic;
lambda1 = S1k\l1;
lambda1theo_Pih = Vh.Pi_h(lambda_theo);
t1 = toc(t1);
show(lambda1theo_Pih);
hold on
show(lambda1);


err0 = sqrt((lambda0 - lambda_theo)|conj((lambda0 - lambda_theo)))/sqrt(lambda0|conj(lambda0))
err1 = sqrt(lambda1 - lambda_theo|conj(lambda1 - lambda_theo))/sqrt(lambda1|conj(lambda1))



























































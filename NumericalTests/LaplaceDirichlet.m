%% Test case on some shape with edges.


Main;

curve = polygonCurve([0,0;0,1;-1,1;-1,0;-1,-1;0 -1;1 -1; 1 0],true,'left');
%unitSegment();
% l = length(curve);
l = 8;
k = 0;
X = R2toRfunc.X;
incWave = 1./sqrt(((X^2 + 0.001^2)));


N = 9;
meshAdapt = MeshCurve(curve,N,@cos,[-pi,0]);
Vh =  weightedFEspace(meshAdapt,'P1','1/sqrt(1-t^2)',...
    'quadNum',3,'specialQuadSegs',1:meshAdapt.nseg);
M = Vh.Mass.concretePart;
[L,U,P,Q] = lu(M);
invM = @(u)(Q*(U\(L \(P*u))));
Wh =  weightedFEspace(meshAdapt,'P1','sqrt(1-t^2)',3);
dM =  Wh.dMass.concretePart;
omega2 = Wh.Mass.concretePart;
K1 = dM - k^2*(omega2 - M);
K =  K1 - k^2*M;

Op_opt = {'tol',1e-3,'fullMatrix',true};
Sw = singleLayer(Vh,...
    'Op_opt',Op_opt,'k',k);
% 
Swgalerk = Sw.galerkine(Vh,'U');    
T0_scal_phi = Vh.phi'*Vh.W;
T0_star_galerk = T0_scal_phi*T0_scal_phi'/sum(Vh.W);

PrecTref = @(u)(invM(TrefethenSqrt(dM,6,invM(u),M,1.5,2*Vh.ndof^2)) + (1/log(2))^2*invM(T0_star_galerk*invM(u)));
secondMemb = Vh.secondMember(-incWave);


t0 = tic;
[lambda0,FLAG0,RELRES0,ITER0,RESVEC0] = variationalSol(Swgalerk,secondMemb,[],1e-8,N,[]);
t0 = toc(t0);
disp(t0);

t1 = tic;
[lambda1,FLAG1,RELRES1,ITER1,RESVEC1] = variationalSol(Swgalerk,secondMemb,[],1e-8,N,PrecTref);
t1 = toc(t1);
disp(t1);


figure
semilogy(1:length(RESVEC0),RESVEC0/norm(secondMemb.concretePart),'-o');
hold on
semilogy(1:length(RESVEC1),RESVEC1/norm(PrecTref(secondMemb.concretePart)),'--x');


xlabel('Iteration number')
ylabel('Residual error')

legend({'Without preconditioner','With preconditioner'});
legend boxoff

set(gca,'FontSize',15)

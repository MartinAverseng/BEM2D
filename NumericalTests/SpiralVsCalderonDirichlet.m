%% Test case on some shape with edges.


Main;
nn = 50;
curve = unitSegment();
l = length(curve);
k = nn*pi/l;
[curve,incWave] = unitSegment(k);


plot(curve);
N = fix(10*k);
meshAdapt = MeshCurve(curve,N,@cos,[-pi,0]);    
Vh =  weightedFEspace(meshAdapt,'P1','1/sqrt(1-t^2)',...
    'quadNum',3,'specialQuadSegs',1:meshAdapt.nseg);
M = Vh.Mass.concretePart;
[L,U,P,Q] = lu(M);
invM = @(u)(Q*(U\(L \(P*u))));
Wh =  weightedFEspace(meshAdapt,'P1','sqrt(1-t^2)',3);

dM =  Wh.dMass.concretePart;
omega2 = Wh.Mass.concretePart;

K1 = dM - k^2*(omega2 -M);
Np = 15;
theta = pi/3;
keps = k+1i*0.025*k^(1/3);
sqrtDarbasK1 = @(x)(padePrecondDarbas(x,Np,theta,keps,M,K1));

Op_opt = {'tol',1e-3,'a_factor',20};
Sw = singleLayer(Vh,...
    'Op_opt',Op_opt,'correcMethod','constantTerm','k',k);

Swgalerk = Sw.galerkine(Vh,'U');

Nw = hyperSingular_w(Vh,'k',k,'Op_opt',Op_opt);
Nwgalerk = Nw.galerkine;



PrecCald = @(u)(omega2\(Nwgalerk*invM(u)));


PrecDarbas = @(u)(invM(sqrtDarbasK1(invM(u))));
clear M L U Q P dM K K1 omega2 

secondMemb = Vh.secondMember(-incWave);

t0 = tic;
[lambda0,FLAG0,RELRES0,ITER0,RESVEC0] = variationalSol(Swgalerk,secondMemb,[],1e-8,N,PrecCald);
t0 = toc(t0);
disp(t0)

t1 = tic;
[lambda1,FLAG1,RELRES1,ITER1,RESVEC1] = variationalSol(Swgalerk,secondMemb,[],1e-8,N,PrecDarbas);
t1 = toc(t1);
disp(t1);


figure
semilogy(1:length(RESVEC0),RESVEC0/norm(PrecCald(secondMemb.concretePart)),'-o');
hold on
semilogy(1:length(RESVEC1),RESVEC1/norm(PrecDarbas(secondMemb.concretePart)),'--x');


xlabel('Iteration number')
ylabel('Residual error')

legend({'Calderon preconditioner','Square root preconditioner'});
legend boxoff


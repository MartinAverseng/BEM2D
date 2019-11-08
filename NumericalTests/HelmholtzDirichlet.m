%% Test case on some shape with edges.


Main;
nn = 10;
curve = spirale;
l = length(curve);
k = nn*pi/l;
[curve,incWave] = spirale(k);


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

Op_opt = {'tol',1e-4,'a_factor',5};
Sw = singleLayer(Vh,...
    'Op_opt',Op_opt,'correcMethod','constantTerm','k',k);

Swgalerk = Sw.galerkine(Vh,'U');

PrecDarbas = @(u)(invM(sqrtDarbasK1(invM(u))));
T0_scal_phi = Vh.phi'*Vh.W;
T0_star_galerk = T0_scal_phi*T0_scal_phi'/sum(Vh.W);

PrecTref = @(u)(invM(TrefethenSqrt(dM,6,invM(u),M,1.5,2*Vh.ndof^2)) + (1/log(2))^2*invM(T0_star_galerk*invM(u)));


% clear M L U Q P dM K K1 omega2 

secondMemb = Vh.secondMember(-incWave);

t0 = tic;
[lambda0,FLAG0,RELRES0,ITER0,RESVEC0] = variationalSol(Swgalerk,secondMemb,[],1e-8,size(Swgalerk,1));
t0 = toc(t0);
disp(t0)

t1 = tic;
[lambda1,FLAG1,RELRES1,ITER1,RESVEC1] = variationalSol(Swgalerk,secondMemb,[],1e-8,size(Swgalerk,1),PrecDarbas);
t1 = toc(t1);
disp(t1);

% t2 = tic;
% [lambda2,FLAG2,RELRES2,ITER2,RESVEC2] = variationalSol(Swgalerk,secondMemb,[],1e-8,50,PrecTref);
% t2 = toc(t2);
% disp(t2);


figure
semilogy(1:length(RESVEC0),RESVEC0/norm(secondMemb.concretePart),'-o');
hold on
semilogy(1:length(RESVEC1),RESVEC1/norm(PrecDarbas(secondMemb.concretePart)),'--x');
%semilogy(1:length(RESVEC2),RESVEC2/norm(PrecTref(secondMemb.concretePart)),'--x');


xlabel('Iteration number')
ylabel('Residual error')

legend({'Without preconditioner','With preconditioner',});
legend boxoff


set(gca,'FontSize',15)

% 

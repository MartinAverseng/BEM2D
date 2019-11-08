%% Test case on the L-shaped domain.

Main;
curve = spirale();

tTot = tic;
N = 400;
meshAdapt = MeshCurve(curve,N,@cos,[-pi,0]);
Vh =  weightedFEspace(meshAdapt,'P1','1/sqrt(1-t^2)',...
    'quadNum',5,'specialQuadSegs',1:meshAdapt.nseg);
t1 = tic;
Op_assemblingOptions = {'tol',1e-4,'a_factor',8};
k = 21;
Sw = singleLayer(Vh,...
    'Op_opt',Op_assemblingOptions,'correcMethod','constantTerm','k',k);
Swgalerk = Sw.galerkine(Vh,'U');
theta_inc = pi;
X =  R2toRfunc(@(Z)(Z(:,1)));
Y =  R2toRfunc(@(Z)(Z(:,2)));
planeWave = exp(1i*k*(X*cos(theta_inc) + Y*sin(theta_inc)));
l = Vh.secondMember(-planeWave);
[lambda2,FLAG2,RELRES2,ITER2,RESVEC2] = variationalSol(Swgalerk,l,[],1e-7,N,Swgalerk.concretePart);
figure
%semilogy(1:length(RESVEC1),RESVEC1,'-o');
%hold on
semilogy(1:length(RESVEC2),RESVEC2,'-o');
drawnow;


%% Compute the radiating solution.

x1 = linspace(-4,4,500);
x2 = linspace(-4,4,500);

[X1,X2] = meshgrid(x1,x2);
Sw.set_X([X1(:) X2(:)]);
valsDiffr = Sw*lambda2;
valsInc = planeWave([X1(:) X2(:)]);
valsDiffr = reshape(valsDiffr,size(X1,1),size(X1,2));
valsInc = reshape(valsInc,size(X1,1),size(X1,2));



%% Animate the wave

amplitude = valsInc+valsDiffr;

figure
imagesc(x1,x2,20*log(abs(amplitude)));
axis xy;
axis equal

tTot = toc(tTot);

figure
animateWave(x1,x2,k,amplitude)



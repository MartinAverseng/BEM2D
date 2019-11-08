%% Test case on some shape with edges.


Main;
nn = 20;
curve = unitSegment;
l = length(curve);
k = nn*pi/l;
[curve,incWave] = unitSegment(k);

plot(curve);
N = fix(10*k);
meshAdapt = MeshCurve(curve,N,@cos,[-pi,0]);    
plot(meshAdapt);
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
keps = k+1i*0.04*k^(1/3);
sqrtDarbasK1 = @(x)(padePrecondDarbas(x,Np,theta,keps,M,K1));

Op_opt = {'tol',1e-3,'a_factor',5};
Sw = singleLayer(Vh,...
    'Op_opt',Op_opt,'correcMethod','constantTerm','k',k);

Swgalerk = Sw.galerkine(Vh,'U');

PrecDarbas = @(u)(invM(sqrtDarbasK1(invM(u))));
% dM12_0 = M*sqrtm(full(M\(-k^2*omega2 + dM)));
% PrecDarbas = @(u)(invM(dM12_0*invM(u)));
dM12_1 = M*sqrtm(full(M\(-k^2*M + dM)));
dM12_2 = M*sqrtm(full(M\(0.2*M + dM)));
Prec1fool = @(u)(invM(dM12_1*invM(u)));
Prec2fool = @(u)(invM(dM12_2*invM(u)));
clear M L U Q P dM K K1 omega2 

secondMemb = Vh.secondMember(-incWave);

t0 = tic;
[lambda0,FLAG0,RELRES0,ITER0,RESVEC0] = variationalSol(Swgalerk,secondMemb,[],1e-8,200);
t0 = toc(t0);
disp(t0)

t1 = tic;
[lambda1,FLAG1,RELRES1,ITER1,RESVEC1] = variationalSol(Swgalerk,secondMemb,[],1e-8,200,PrecDarbas);
t1 = toc(t1);
disp(t1);

t2 = tic;
[lambda2,FLAG2,RELRES2,ITER2,RESVEC2] = variationalSol(Swgalerk,secondMemb,[],1e-8,200,Prec1fool);
t2 = toc(t2);
disp(t2);

t3 = tic;
[lambda3,FLAG3,RELRES3,ITER3,RESVEC3] = variationalSol(Swgalerk,secondMemb,[],1e-8,200,Prec2fool);
t3 = toc(t3);
disp(t3);

figure
semilogy(1:length(RESVEC0),RESVEC0/norm(secondMemb.concretePart),'-o');
hold on
semilogy(1:length(RESVEC3),RESVEC3/norm(Prec2fool(secondMemb.concretePart)),'--x');
semilogy(1:length(RESVEC2),RESVEC2/norm(Prec1fool(secondMemb.concretePart)),'--x');
semilogy(1:length(RESVEC1),RESVEC1/norm(PrecDarbas(secondMemb.concretePart)),'--x');


xlabel('Iteration number')
ylabel('Residual error')

legend({'Pas de preconditioneur','$\sqrt{-(\omega\partial_x)^2 + 0.2I_d}$','$\sqrt{-(\omega\partial_x)^2 -k^2I_d}$','$\sqrt{-(\omega \partial_x)^2 - k^2 \omega^2}$'},'Interpreter','latex');
legend boxoff


figure
x1 = -2:0.01:2;
x2 = -2:0.01:2;
[Xfig,Yfig] = meshgrid(x1,x2);
Z = [Xfig(:) Yfig(:)];
RadiateS = singleLayer(Vh,'X',Z,'k',k,'Op_opt',{'a_factor',4,'tol',1e-3});
uD = Vh.secondMember(-incWave);
uS = RadiateS*lambda1;
uinc = incWave(Z);
uinc = reshape(uinc,size(Xfig,1),size(Xfig,2));
uS = reshape(uS,size(Xfig,1),size(Xfig,2));
amplitude= (uinc + uS);
amplitude = reshape(amplitude,size(Xfig,1),size(Xfig,2));
figure
imagesc(x1,x2,20*log(abs(amplitude)));
hold on
plot(curve)
axis xy;
axis equal
colormap gray
caxis([-5 6]);
hold on;
plot(curve);
axis tight



t = 0;
dt = 0.005;
figure 
plot(curve,'k');
% id = 1/a^2*Xfig.^2 + + 1/b^2*Yfig.^2 <=1;
while true
    t = t + dt;
    toPlot = 0;
    for j = 1:length(k)
        amplitude = amplitude*exp(-1i*k*dt);
        toPlot = toPlot + amplitude;
    end
%     toPlot(id) = nan;
    imagesc(x1,x2,real(toPlot));
%     plot(curve,'k');
    axis equal
    axis xy
    caxis([-3,3])
    drawnow;
end

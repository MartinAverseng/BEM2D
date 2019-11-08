clear all;
close all;

curve = archSpirale();
l = length(curve);
nn = 50;
k0 = nn*pi/l;
N = fix(50*k0);
mesh = MeshCurve(curve,N,@cos,[-pi,0]);

Vh = weightedFEspace(mesh,'P1','1/sqrt(1-t^2)',3);
plot(mesh);

Wh =  weightedFEspace(mesh,'P1','sqrt(1-t^2)',3);

dM =  Wh.dMass.concretePart;
omega2 = Wh.Mass.concretePart;

x0 = -5;
k = linspace(k0-2,k0+2,40);
dk = k(2) - k(1);
% dk = 1;
% k = [k,k0];
x1 = -2:0.01:2;
x2 = -2:0.01:2;
[Xfig,Yfig] = meshgrid(x1,x2);
Z = [Xfig(:) Yfig(:)];

M = Vh.Mass.concretePart;
[L,U,P,Q] = lu(M);
invM = @(u)(M\u);%@(u)(Q*(U\(L \(P*u))));
% dM =  Vh.dMass.concretePart;
% Np = 15;
% theta = pi/3;
% keps = k+1i*0.01*k^(1/3);
% sqrtDarbasK1 = @(x)(padePrecondDarbas(x,Np,theta,keps,M,dM));
% clear dM;
% clear omega2;
% u = randn(N+1,1);
% plot((dM - k^2*omega2)*u);
% hold on
% plot(real(sqrtDarbasK1(M\sqrtDarbasK1(u))));

theta_inc = 0;

for j = 1:length(k)
    ak(j) = exp(-(k(j) - k0)^2/2)*exp(-1i*k(j)*x0)*dk;
    Sw = singleLayer(Vh,'k',k(j),'Op_opt',{'fullMatrix',true});
    nZ2 = Z(:,1).^2 + Z(:,2).^2;
%     idx = and(and(nZ2>0.6.^2,nZ2<1.4^2),Z(:,1)>-0.6);
    RadiateSw = singleLayer(Vh,'X',Z,'k',k(j),'Op_opt',{'a_factor',10,'tol',1e-3});
%     RadiateSw2 = singleLayer(Vh,'X',Z(~idx,:),'k',k(j),'Op_opt',{'a_factor',60,'tol',1e-3});
    Swh = Sw.galerkine(Vh,'U');
    
    K1 = dM - k(j)^2*(omega2 -M);
    Np = 15;
    theta = pi/3;
    keps = k(j)+1i*0.025*k(j)^(1/3);
    sqrtDarbasK1 = @(x)(padePrecondDarbas(x,Np,theta,keps,M,K1));
    
    PrecDarbas = @(u)(invM(sqrtDarbasK1(invM(u))));
    
    PW = R2toRfunc(@(Z)(exp(1i*((k(j)))*(Z(:,1)*cos(theta_inc) + Z(:,2)*sin(theta_inc)))));
    %     PW = R2toRfunc(@(Z)(1./(abs(Z(:,1)-a)+ eps)));
    uD = Vh.secondMember(-PW);
    [lambda,flag,relres,iter,resvec] = variationalSol(Swh,uD,[],1e-10,25,Swh.concretePart);
    %     [lambda2,flag2,relres2,iter2,resvec2] = variationalSol(Sh,uD,[],1e-10,size(Sh,1),Prec);
%     lambda = Swh.concretePart\uD.concretePart;
%     lambda = FE_func(Vh,lambda)
%     tmp = 0*Z(:,1);
%     tmp(idx) = RadiateSw1*lambda;
%     tmp(~idx) = RadiateSw2*lambda;
    uS{j} = RadiateSw*lambda;
    uinc{j} = PW(Z);
    uinc{j} = ak(j)*reshape(uinc{j},size(Xfig,1),size(Xfig,2));
    uS{j} = ak(j)*reshape(uS{j},size(Xfig,1),size(Xfig,2));
    amplitude{j} = (uinc{j} + uS{j});
    amplitude{j} = reshape(amplitude{j},size(Xfig,1),size(Xfig,2));
%     figure
%     imagesc(x1,x2,20*log(abs(amplitude{j})));
%     hold on
%     plot(curve)
%     axis xy;
%     axis equal
    fprintf('*');
end
disp(' ');

% figure
% semilogy(1:iter(end)+1,resvec/norm(uD.concretePart));
% hold on
% semilogy(1:iter2(end)+1,resvec2/norm(Prec(uD.concretePart)));
% drawnow;

clear Sw;
clear Swh;
clear uS
clear uinc
clear radiateSw;
clear Z;
clear Xfig;
clear Yfig;
close all;


t = 0;
dt = 0.2;
% id = 1/a^2*Xfig.^2 + + 1/b^2*Yfig.^2 <=1;
figure
toPlot = 0;
for j = 1:length(k)
    amplitude{j} = amplitude{j}*exp(-1i*k(j)*dt);
    toPlot = toPlot + amplitude{j};
end
axis equal
axis xy
drawnow;

while true
    t = t + dt;
    toPlot = 0;
    for j = 1:length(k)
        amplitude{j} = amplitude{j}*exp(-1i*k(j)*dt);
        toPlot = toPlot + amplitude{j};
    end
    %     toPlot(id) = nan;
    imagesc(x1,x2,(20*log(abs(toPlot))));
    %     plot(curve,'k');
    caxis([-50,100])
    axis equal
    axis xy
    drawnow;
    colormap gray
    
end



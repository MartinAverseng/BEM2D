clear all;
close all;

a = 1;
b = 2;
LshapedVertices = [0 0; 1 0; 1 1; 1/2 1; 1/2 1/2; 0 1/2];
curve = polygonCurve(LshapedVertices,true,'left');
curve = ellipse(a,b);
% curve = boomerang();
N = 100;
mesh = MeshCurve(curve,N);

Vh = FEspace(mesh,'P1',3);
plot(mesh);

k0 = 5;
kR = k0*curve.length;
x0 = -5;
k = linspace(k0-1,k0+1, 10);
dk  = k(2) - k(1);
Dk = 1/(k(end) - k(1));
% k = [k,k0];
figure
x1 = -3:0.05:3;
x2 = -3:0.05:3;
[Xfig,Yfig] = meshgrid(x1,x2);
Z = [Xfig(:) Yfig(:)];

M = Vh.Mass.concretePart;
[L,U,P,Q] = lu(M);
invM = @(u)(M\u);%@(u)(Q*(U\(L \(P*u))));
dM =  Vh.dMass.concretePart;


thetai = -pi/6;
lambda0 = FE_func(Vh,0);
for j = 1:length(k)
    ak(j) = exp(-2*(k(j) - k0)^2/2)*exp(-1i*k(j)*x0)*dk/Dk;
    tic; S = singleLayer(Vh,'k',k(j),'Op_opt',{'fullMatrix',true,'tol',1e-5});toc;
%     tic; H = hyperSingular(Vh,'k',k(j),'Op_opt',{'fullMatrix',true,'tol',1e-5});toc;
%     tic; D = doubleLayer(Vh,'k',k(j),'Op_opt',{'fullMatrix',true}); toc;
    Sh = S.galerkine(Vh,'U');
%     Dh = D.galerkine(Vh,'U');
%     Hh = H.galerkine;
%     norm(Sh.concretePart*M^(-1)*Hh.concretePart - M/4 - Dh.concretePart*M^(-1)*Dh.concretePart)
    %    Dh = D.galerkine(Vh,'U');
    keps = k(j) + 1i*0.05;
%     Y = 1/2*(dM-k(j)^2*M);
    Np = 40;
    theta = pi/3;
    keps = k(j)+1i*0.01*k(j)^(1/3);
    sqrtDarbasK1 = @(x)(padePrecondDarbas(x,Np,theta,keps,M,dM));
    PrecDarbas = @(u)(M\(sqrtDarbasK1(M\u)));
% %     clear dM;
%     u = randn(N,1);
%     plot((dM - k(j)^2*M)*u);
%     hold on
%     plot(real(sqrtDarbasK1(M\sqrtDarbasK1(u))));
%     [P,DD] = eig(full(M\Y));
    %     i = find(abs(diag(DD)) < 1e-6);
%     %     D(i,i) = 1;
%     sqrtdM = M*P*sqrt(DD)*P^(-1);
    %     Mat = Sh.concretePart*M^(-1)*sqrtdM + M/2 + Dh.concretePart;
    %     bili2 = BilinearForm(Vh,Vh,Mat);
    %     disp(cond(M\Sh.concretePart));
    %     disp(cond(M\(Mat*Mat')));
    PW = R2toRfunc(@(Z)(exp(1i*((k(j)))*(cos(thetai)*Z(:,1) + sin(thetai)*Z(:,2)))));
    %     PW = R2toRfunc(@(Z)(1./(abs(Z(:,1)-a)+ eps)));
    uD = Vh.secondMember(-PW);
    mu = Vh.Pi_h(-PW);
    %     [lambda,flag,relres,iter,resvec] = gmres(Sh.concretePart,uD.concretePart,[],1e-10,32);
    %     [lambda2,flag2,relres2,iter2,resvec2] = gmres(Mat'*Mat,Mat'*uD.concretePart,[],1e-10,32);
    %     figure
    %     semilogy(1:iter(end)+1,resvec/norm(uD.concretePart),'-o');
    %     hold on
    %     semilogy(1:iter2(end)+1,resvec2/norm(Mat'*uD.concretePart),'-x');
    %     xlabel("Nombre d'itérations")
    %     ylabel("Résidu");
    %     title('Convergence de GMRES')
    %     legend({'Système original','Système préconditionné'})
    %     legend('Location','Southwest')
    %     drawnow;
    %     lambda = variationalSol(Sh,(1/2)*uD - Dh*mu,[],1e-10,[],Sh.concretePart);
%     PrecCald = @(u)(M\(Hh*(M\u)));

    [lambda,~,~,iter1,resvec1] = variationalSol(Sh,uD,[],1e-10,size(Sh,1));
    [lambda2,~,~,iter2,resvec2] = variationalSol(Sh,uD,[],1e-10,size(Sh,1),PrecDarbas);
    figure;
    semilogy(1:length(resvec1),resvec1/norm(uD.concretePart),'--x');
    hold on
    semilogy(1:length(resvec2),resvec2/norm(PrecDarbas(uD.concretePart)),'--x');
    
    % lambda = FE_func(Vh,lambda);
    RadiateS = singleLayer(Vh,'X',Z,'k',k(j),'Op_opt',{'a_factor',8,'tol',1e-6});
    %     RadiateD = doubleLayer(Vh,'X',Z,'k',k(j),'Op_opt',{'fullMatrix',true,'tol',1e-5});
    
    uS{j} = RadiateS*(lambda); %+ RadiateD*(mu);
    uinc{j} = PW(Z);
    uinc{j} = ak(j)*reshape(uinc{j},size(Xfig,1),size(Xfig,2));
    uS{j} = ak(j)*reshape(uS{j},size(Xfig,1),size(Xfig,2));
    amplitude{j} = (uinc{j} + uS{j});
    amplitude{j} = reshape(amplitude{j},size(Xfig,1),size(Xfig,2));
    %     figure
    %     imagesc(x1,x2,real((uS{j})));
    %     hold on
    %     plot(curve)
    %     axis xy;
    %     axis equal
end

close all;
clear RadiateS;
clear RadiateD;
clear S;
clear Sh;
clear Dh;
clear S;
clear D;
clear uS;

t = 0;
dt = 0.05;
% id = 1/a^2*Xfig.^2 + + 1/b^2*Yfig.^2 <=1;
figure
while true
    t = t + dt;
    toPlot = 0;
    for j = 1:length(k)
        amplitude{j} = amplitude{j}*exp(-1i*k(j)*dt);
        toPlot = toPlot + amplitude{j};
    end
    %     toPlot(id) = nan;
    imagesc(x1,x2,real(toPlot));
    %     plot(curve,'k');
    axis equal
    axis xy
    caxis([-2,2])
    drawnow;
end

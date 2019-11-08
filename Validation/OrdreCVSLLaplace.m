%% Ordre de CV Laplace sur un disque de rayon 1/2 (pour assurer que S est positif).
clear all;
close all;


Ns = [20, 50, 100, 200, 500, 1000, 5000];

a = 1/4;
b = 1/4;
curve = ellipse(a,b);
X = R2toRfunc.X; Y = R2toRfunc.Y;
phi = atan2(Y,X);
m =5;
u0 = cos(m*phi);

lambda = 4*(2*m)*u0;
k = 0;
errH_12 = zeros(length(Ns),1);
errL2 = zeros(length(Ns),1);
for i = 1:length(Ns)
    N = Ns(i)
    
    mesh = MeshCurve(curve,N);
    
    Vh = FEspace(mesh,'P0',3);
    M = Vh.Mass.concretePart;
    
    S = singleLayer(Vh,'k',k,'Op_opt',{'fullMatrix',true});
    Swgalerk = S.galerkine(Vh,'U');
    
    u0_h = Vh.secondMember(u0);
    
    lambda_h = variationalSol(Swgalerk,u0_h,[],1e-10,[],Swgalerk.concretePart);
    eh = Vh.Pi_h(lambda) - lambda_h;
    errH_12P0(i) = sqrt(Swgalerk*eh | eh);
    errL2P0(i) = sqrt(real((lambda - lambda_h)|(lambda - lambda_h)));
end

for i = 1:length(Ns)
    N = Ns(i)
    
    mesh = MeshCurve(curve,N);
    
    Vh = FEspace(mesh,'P1',3);
    M = Vh.Mass.concretePart;
    
    S = singleLayer(Vh,'k',k,'Op_opt',{'fullMatrix',true});
    Swgalerk = S.galerkine(Vh,'U');
    
    u0_h = Vh.secondMember(u0);
    
    lambda_h = variationalSol(Swgalerk,u0_h,[],1e-10,[],Swgalerk.concretePart);
    eh = lambda - lambda_h;
    
    errL2P1(i) = sqrt(real((lambda - lambda_h)|(lambda - lambda_h)));
end

hs = pi*(a + b)./Ns;
loglog(hs,errL2P0,'DisplayName','$P^0$');
hold on
loglog(hs,hs/hs(1)*errL2P0(1),'k--','HandleVisibility','off');
loglog(hs,errL2P1,'DisplayName','$P^1$');
hold on
loglog(hs,hs.^(2)/hs(1)^(2)*errL2P1(1),'k--','DisplayName','Théorie');


title('Ordres de convergence numériques')
xlabel('$h$','Interpreter','latex');
ylabel('$||\lambda - \lambda_h ||_{L^2(\Gamma)}$','Interpreter','latex')
set(gca,'FontSize',15);
L = legend;
set(L,'Interpreter','latex');
legend show
legend('location','Southeast')

xlim([hs(end),hs(1)])
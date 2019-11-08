% Validate the hypersingular potential on the circle.

Main;
N = 1000;
q = 3;

curve = circle(1,[0,0]);
mesh = MeshCurve(curve,N);

Vh = FEspace(mesh,'P1');
Wh = FEspace(mesh,'P2');

N1 = Vh.hyperSingularOperator('ZeroMeanModif',0.001);
N2 = Wh.hyperSingularOperator('ZeroMeanModif',0.001);

u0_func = @(Z)(Z(:,1).^6 - 15*Z(:,1).^4.*Z(:,2).^2 + 15*Z(:,1).^2.*Z(:,2).^4 - Z(:,2).^6);
u0 = R2toRfunc(u0_func);
lambda_theo = 1/3*u0;
lambda_theo_1 = Vh.Pi_h(lambda_theo);
lambda_theo_2 = Wh.Pi_h(lambda_theo);

l1 = Vh.secondMember(u0);
l2 = Wh.secondMember(u0);
test = (1|lambda_theo_1);

t1 = tic;
lambda1 = variationalSol(N1,l1);
t1 = toc(t1);
t2 = tic;
lambda2 = variationalSol(N2,l2);
t2 = toc(t2);

% showLog(lambda1 - lambda_theo);
% hold on
% showLog(lambda2 - lambda_theo);

err1 = sqrt(lambda1-lambda_theo|lambda1-lambda_theo);
disp('P1 L^2 error:');
disp(err1);
err2 = sqrt(lambda2-lambda_theo|lambda2-lambda_theo);
disp('P2 L^2 error:');
disp(err2);

err12 = sqrt(u0|lambda_theo-lambda1);
disp('P1 H^{1/2} error:')
disp(err12);
err22 = sqrt(u0|lambda_theo-lambda2);
disp('P2 H^{1/2} error:')
disp(err22);

test1 = sqrt(N1*(lambda1-lambda_theo_1)|(lambda1-lambda_theo_1));
test2 = sqrt(N2*(lambda2-lambda_theo_2)|(lambda2-lambda_theo_2));




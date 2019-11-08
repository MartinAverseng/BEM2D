% Validate double layer
% Check Gauss identity
clear all;
close all;
curve = ellipse(1,1);
plot(curve)
N = 100;
mesh = MeshCurve(curve,N);
plot(mesh);
Vh = FEspace(mesh,'P1',3);
% 
x1 = -10:0.1:10;
x2 = -10:0.1:10;
[X,Y] = meshgrid(x1,x2);
Z = [X(:),Y(:)];
S = singleLayer(Vh,'Xdata',Z,'Op_opt',{'fullMatrix',true});
D = doubleLayer(Vh,'Xdata',Z,'Op_opt',{'fullMatrix',true});

figure
lambda = FE_func(Vh,1);
vals = D*lambda;
vals = reshape(vals,size(X,1),size(X,2));
figure
imagesc(x1,x2,((vals)));
axis xy;
axis equal

figure
lambda = FE_func(Vh,1);
vals = S*lambda;
vals = reshape(vals,size(X,1),size(X,2));
figure
imagesc(x1,x2,((vals)));
axis xy;
axis equal

% Check spectrum 
D = doubleLayer(Vh,'Op_opt',{'fullMatrix',true});
Dh = D.galerkine(Vh,'U');
M = Vh.Mass;
[P,D] = eig(M\Dh.concretePart);
disp(diag(D));


% % Double layer at frequency k :
% k = 15;
% D = doubleLayer(Vh,'k',k,'Xdata',Z,'Op_opt',{'fullMatrix',true});
% lambda = FE_func(Vh,1);
% vals = D*lambda;
% vals = reshape(vals,size(X,1),size(X,2));
% figure
% imagesc(x1,x2,real(vals));
% axis xy;
% axis equal
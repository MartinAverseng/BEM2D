% 
clear all;
close all;
curve = ellipse(1,1);
plot(curve)
N = 30;
mesh = MeshCurve(curve,N);
plot(mesh);
Vh = FEspace(mesh,'P1',3);
M = Vh.Mass;
% 
x1 = -10:0.1:10;
x2 = -10:0.1:10;
[X,Y] = meshgrid(x1,x2);
Z = [X(:),Y(:)];
H = hyperSingular(Vh,'k',2,'Op_opt',{'fullMatrix',true});
Hh = H.galerkine;
[P,D] = eig(M\Hh.concretePart);
disp(sort(diag(D)));
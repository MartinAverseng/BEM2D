function [N,AB,XB,XA,l,u,a,b,d] = parameters_singInt(X,A,B)

N = normalVector(A,B);
AB = B-A;
XB = B - X;
XA = A - X;
l = sqrt(cWise_dot(AB,AB));
u = AB./[l l];
a = cWise_dot(XA,u);
b = cWise_dot(XB,u);
d = abs(cWise_dot(N,XA));



function[Z] = cWise_dot(X,Y)
    Z = X(:,1).*Y(:,1) + X(:,2).*Y(:,2);
end

end


function [ out ] = J_1( X,A,B,N )
% computes int_{a}^b (x-y).n/|x-y|^2 dy where n is the normal vector of the
% segment [A,B].

AB = B-A;
XB = B - X;
XA = A - X;
l = sqrt(cWise_dot(AB,AB));
u = AB./[l l];
a = cWise_dot(XA,u);
b = cWise_dot(XB,u);
d = cWise_dot(N,-XA);

F = @(x)(primitive(x));
out = F(b)-F(a);

    function[Z] = cWise_dot(X,Y)
        Z = X(:,1).*Y(:,1) + X(:,2).*Y(:,2);
    end
    function[out] = primitive(x)
        out = 0*x;
        ind = abs(d)>1e-10;
        out(ind) = 1/2*d(ind).*log(d(ind).^2+x(ind).^2);
        out(~ind) = 0;
    end

end


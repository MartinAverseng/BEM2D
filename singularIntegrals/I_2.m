function [ out ] = I_2(X,A,B)

% This function computes the value of 
% \int_[A,B] \ln|x-y|y^2 dy
% where X is any point in the plane, [A,B] is any segment, and the optional
% argument N is a unit normal vector to the segment [A,B] used for the
% computation. 
    
[~,~,~,~,~,~,a,b,d] = parameters_singInt(X,A,B);

F=@(x)((1/6)*x.^3.*log(d.^2+x.^2)-(1/9)*x.^3+(1/3)*d.^2.*x-(1/3)*d.^3.*atan(x./d));
out = F(b) - F(a);


end


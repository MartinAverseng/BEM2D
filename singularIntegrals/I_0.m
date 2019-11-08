function [ out ] = I_0(X,A,B)

% This function computes the value of 
% \int_[A,B] \ln|x-y|dy
% where X is any point in the plane, [A,B] is any segment, and the optional
% argument N is a unit normal vector to the segment [A,B] used for the
% computation. 

[~,~,~,~,~,~,a,b,d] = parameters_singInt(X,A,B);
F = @(x)(1/2*x.*log(d.^2 + x.^2)-x + d.*atan(x./d));
out = F(b) - F(a);

end

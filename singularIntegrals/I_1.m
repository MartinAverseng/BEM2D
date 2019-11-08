function [ out ] = I_1(X,A,B)

% This function computes the value of 
% \int_[A,B] \ln|x-y| |y-A| dy
% where X is any point in the plane, [A,B] is any segment, and the optional
% argument N is a unit normal vector to the segment [A,B] used for the
% computation. 
    
[~,~,~,~,~,~,a,b,d] = parameters_singInt(X,A,B);
F = @(x)(1/4*(x.^2 + d.^2).*(log(x.^2 + d.^2)-1));
out = F(b) - F(a);


end

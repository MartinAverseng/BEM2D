function [ N ] = normalVector(A,B)
% returns a normal unit norm vector of the segment [A,B]
% A and and B are n x 2 arrays, and N is returned as a n x 2 array

AB = B-A;
l = sqrt(cWise_dot(AB,AB));
T = [AB(:,1)./l AB(:,2)./l];
Rot = [0 1; -1 0];
N = T*Rot;



end


function[Z] = cWise_dot(X,Y)
    Z = X(:,1).*Y(:,1) + X(:,2).*Y(:,2);
end

function[B] = squarePol(A)
% A is a matrix with each column representing the coefficients
% of a polynom.
r = size(A,1);
M = size(A,2);
l = 2*r-1;
B = zeros(l,M);
for m = 1:M
    B(:,m) = conv(A(:,m),A(:,m));
end
end
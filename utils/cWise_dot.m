function[Z] = cWise_dot(X,Y)
    Z = X(:,1).*Y(:,1) + X(:,2).*Y(:,2);
end
function [ I ] = isClose(X,Y,rmin)

[I,~] = rangesearch(X,Y,rmin);

end


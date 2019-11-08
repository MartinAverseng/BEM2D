function [ out ] = xlog(x)

out = x.*log(abs(x));
out(isnan(out)) = 0;

%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

end


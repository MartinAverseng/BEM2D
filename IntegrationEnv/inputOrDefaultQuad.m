function [ q ] = inputOrDefaultQuad(q)

if ~exist('q','var')||isempty(q)
    q = 3;
end
if isa(q,'double')
    q = NumericalQuadrature(q);
else
    assert(isa(q,'NumericalQuadrature'));
end


end


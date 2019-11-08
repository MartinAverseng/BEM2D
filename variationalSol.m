function [fe_func,flag,relres,iter,resvec,t] = variationalSol(b,l,varargin)
% Solves the variational problem b(u,v) = l(v)
% for all v in the space of functions on which both b and l are defined.
% returns a FE_func object containing the coordinates of this solution

assert(isa(b,'BilinearForm'));
assert(isa(l,'LinearForm'));
% assert(isequal(b.Wh,l.feSpace));


B = AbstractMatrix(b);
t = tic;
[x,flag,relres,iter,resvec] = B.mldivide(l.concretePart,varargin{:});
t = toc(t);
fe_func = FE_func(b.Vh,x);
%fprintf('\nVariational equation : gmres returned a solution in %s iteration\n',num2str(length(resvec)));


end


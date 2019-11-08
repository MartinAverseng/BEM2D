classdef BilinearForm < AbstractMatrix
    % a bilinear form is a function defined on a Geometry and that return to
    % every pair of FE_func object (u,v) defined on the same geometry a real value
    % approximately equal to \int_{mesh}\int_{mesh} G(x,y)u(x)v(y)dy
    
    properties
        Vh@FEspace;
        Wh@FEspace
    end
    
    methods
        function[this] = BilinearForm(VVh,WWh,varargin)
            this@AbstractMatrix(varargin{:});
            assert(isequal(size(this),[VVh.ndof,WWh.ndof]));
            this.Vh = VVh;
            this.Wh = WWh;
        end
        function[C] = mtimes(A,B)
            if isa(A,'BilinearForm')
                if isa(B,'FE_func')
                    assert(isequal(B.feSpace, A.Wh));
                    C = LinearForm(A.Vh,A*(B.v));
                elseif isa(B,'double')
                    C = mtimes@AbstractMatrix(A,B);
                else
                    error(...
                        ['wrong type of argument, expected ''FE_func'' or scalar ''double'','...
                        'received %s instead'],...
                        class(B));
                end
            else
                % B is then a Bilinear form, so A must be a scalar double.
                if isa(A,'double')
                    C = BilinearForm(B.Vh,B.Wh,A*AbstractMatrix(B));
                else
                    error('Product not possible')
                end
            end
        end
        function[out] = mldivide(this,l)
            if and(isa(this,'BilinearForm'),isa(l,'LinearForm'))
                out = variationalSol(this,l);
            else
                out = mldivide@AbstractMatrix(this,l);
            end
        end
        function[out] = eval(this,uh,vh)
            VVh = this.Vh;
            WWh = this.Wh;
            assert(VVh.contains(uh)); assert(WWh.contains(vh));
            out = sum((this*uh.v).*conj(vh.v));
        end
    end
    
end






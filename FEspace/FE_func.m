classdef FE_func < func
    % Object defining a function f defined on a mesh that is of the form
    % f(x) = \sum(\alpha_i \phi_i(x))
    % where \phi_i is the i-th basis function of the FE space.
    % This is a "pure" FE function, so the handle part is unset. If user
    % tries to set the handle part of a FE_func object, it will produce an
    % error. FE_func can be added with FE_func or func objects, with the
    % following rules :
    % FE_func + FE_func = FE_func
    % FE_func + func = func
    
    properties
    end
    
    methods
        function[this] = FE_func(Vh,vv)
            assert(isa(vv,'double'));
            if isscalar(vv)
                vv = vv*ones(Vh.ndof,1);
            elseif isempty(vv)
                vv = sparse([],[],[],Vh.ndof,1);
            else
                assert(isequal(size(vv),[Vh.ndof 1]));
            end
            this.feSpace = Vh;
            this.fePart = vv;
        end
        function[out] = size(this,dim)
            out = [this.feSpace.ndof,1];
            if nargin == 2
                
                out = out(dim);
                
            end
        end
        function[vv] = v(this)
            vv = this.fePart;
        end
        function[C] = plus(A,B)
            [A_new,B_new] = operationInterpreter(A,B);
            if isequal(class(A_new),'FE_func')
                C = FE_func(commonSpace(A_new,B_new),A_new.fePart + B_new.fePart);
            else
                C = plus@func(A_new,B_new);
            end
        end
        function[A] = conj(A)
            A.fePart = conj(A.fePart);
        end
        function[C] = uminus(A)
            C = -1*A;
        end
        function[C] = times(A,B)
            C = FE_func(commonSpace(A,B),A.fePart.*B.fePart);
        end
        function[C] = rdivide(A,B)
            C = FE_func(commonSpace(A,B),A.fePart./B.fePart);
        end
        function[C] = mtimes(A,B)
            assert(or(isa(A,'double'),isa(B,'double')),sprintf(...
                'One of the arguments should be a scalar double, received istead %s and %s'...
                ,class(A),class(B)));
            if isa(A,'double')
                C = mtimes(B,A);
            else
                assert(and(isa(B,'double'),isscalar(B)));
                C = FE_func(A.feSpace,A.fePart*B);
            end
        end
        function[this] = real(this)
            this.fePart = real(this.v);
        end
        function[this] = imag(this)
            this.fePart = imag(this.v);
        end
        function[this] = modulus(this)
            this.fePart = abs(this.v);
        end
    end
    
end

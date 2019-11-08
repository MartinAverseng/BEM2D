classdef R2toRfunc < func
    % Function R^2 -> R
    properties
    end
    
    methods
        function[this] = R2toRfunc(varargin)
            if nargin == 0
                this.handlePart = @(Z)(sparse(size(Z,1),1)); % Identical 0 func
            end
            if nargin==1
                arg = varargin{1};
                if isa(arg,'function_handle')
                    assert(isequal(size(arg(zeros(2,2))),[2,1]));
                    this.handlePart = arg;
                elseif isa(arg,'double')
                    if isempty(arg)
                        this.handlePart = @(Z)(sparse(size(Z,1),1));
                    else
                        assert(isscalar(arg));
                        this.handlePart = @(Z)(0*Z(:,1) + arg);
                    end
                elseif isa(arg,'R2toRfunc')
                    this = arg;
                end
            end
        end
        function[v] = subsref(this,Z)
            assert(length(Z.subs)==1,'too many input arguments');
            assert(size(Z.subs{1},2)==2,...
                'function must be evaluated on Nx2 args');
            v = this.handlePart(Z.subs{1});
        end
        function[C] = plus(A,B)
            [A_new,B_new] = operationInterpreter(A,B);
            if ~isequal(class(A_new),'R2toRfunc')
                C = plus@func(A_new,B_new);
            else
                C = R2toRfunc(@(Z)(A_new.handlePart(Z) + B_new.handlePart(Z)));
            end
        end
        function[C] = conj(A)
            C = R2toRfunc(@(Z)(conj(A.handlePart(Z))));
        end
        function[C] = uminus(A)
            C = -1*A;
            
        end
        function[C] = mtimes(A,B)
            if and(isa(A,'double'),isscalar(A))
                C = R2toRfunc(@(Z)(A*B.handlePart(Z)));
            elseif and(isa(B,'double'),isscalar(B))
                C = R2toRfunc(@(Z)(A.handlePart(Z)*B));
            else
                assert(isequal(class(A),'R2toRfunc'));
                assert(isequal(class(B),'R2toRfunc'));
                C = R2toRfunc(@(Z)(A.handlePart(Z).*B.handlePart(Z)));
            end
        end
        function[C] = mpower(A,B)
            if and(isa(A,'double'),isscalar(A))
                C = R2toRfunc(@(Z)(A.^(B.handlePart(Z))));
            elseif and(isa(B,'double'),isscalar(B))
                C = R2toRfunc(@(Z)(A.handlePart(Z).^B));
            else
                assert(isequal(class(A),'R2toRfunc'));
                assert(isequal(class(B),'R2toRfunc'));
                C = R2toRfunc(@(Z)(A.handlePart(Z).^(B.handlePart(Z))));
            end
        end
        function[C] = times(A,B)
            [A_new,B_new] = operationInterpreter(A,B);
            if isequal(class(A_new),'R2toRfunc')
                C = mtimes(A,B);
            else
                C = times@func(A_new,B_new);
            end
        end
        function[C] = rdivide(A,B)
            [A_new,B_new] = operationInterpreter(A,B);
            C = R2toRfunc(@(Z)(A_new.handlePart(Z)./B_new.handlePart(Z)));
        end
        function[C] = mrdivide(A,B)
            C = rdivide(A,B);
        end
        function[B] = applyFun(A,arbitraryFun)
            B = R2toRfunc(@(Z)(arbitraryFun(A.handlePart(Z))));
        end
        function[B] = exp(A)
            B = applyFun(A,@exp);
        end
        function[B] = sin(A)
            B = applyFun(A,@sin);
        end
        function[B] = cos(A)
            B = applyFun(A,@cos);
        end
        function[B] = sqrt(A)
            B = applyFun(A,@sqrt);
        end
        function[B] = ln(A)
            B = applyFun(A,@log);
        end
        function[B] = log(A)
            B = applyFun(A,@log);
        end
        function[C] = power(A,n)
            C = R2toRfunc(@(Z)(A.handlePart(Z).^n));
        end
        function[C] = atan2(Y,X)
            C = R2toRfunc(@(Z)(atan2(Y.handlePart(Z),X.handlePart(Z))));
        end
        function[B] = asin(X)
            B = applyFun(X,@asin);
        end
        function[B] = acos(X)
            B = applyFun(X,@acos);
        end
        function[B] = dilog(A)
            B = applyFun(A,@dilog);
        end
        
        function[] = disp(~)
            fprintf('object of class R2toRfunc \n');
        end
    end
    methods (Static)
        function[B] = X
            B = R2toRfunc(@(Z)(Z(:,1)));
        end
        function[B] = Y
            B = R2toRfunc(@(Z)(Z(:,2)));
        end
        function[B] = Tn(n)
            B = R2toRfunc(@(Z)(chebyshevT(n,Z(:,1))));
        end
        function[B] = Un(n)
            B = R2toRfunc(@(Z)(chebyshevU(n,Z(:,1))));
        end
    end
end


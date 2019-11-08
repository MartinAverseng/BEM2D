classdef LinearForm < AbstractMatrix
    % a linear form is a function defined on a FE space, that accepts
    % FE_func arguments defined on the same FE space and returns a scalar.
    
    properties
        feSpace@FEspace;
    end
    
    methods
        function[this] = LinearForm(Vh,values)
            this = this@AbstractMatrix;
            if nargin >0
                this.concretePart = values(:);
                this.feSpace = Vh;
                this.N1 = Vh.ndof;
                this.N2 = 1;
            end
        end
        function[r] = eval(this,u)
            assert(isa(u,'FE_func'));
            errorMessage = 'This linear form is not defined on the same FE space as argument ''u''';
            assert(isequal(u.feSpace,this.feSpace),errorMessage);
            r = sum(this.concretePart.*conj(u.v));
        end
        function[out] = or(this,u)
            out = this.eval(u);
        end
        function[C] = plus(this,l)
            if isequal(class(this),'LinearForm')
                switch class(l)
                    case 'LinearForm'
                        assert(isequal(this.feSpace,l.feSpace),...
                            'The linear forms are not defined on the same FE space')
                        lnew = l;
                    case 'func'
                        lnew = this.feSpace.secondMember(l);
                    case 'FE_func'
                        lnew = this.feSpace.secondMember(l);
                    case 'R2toRfunc'
                        lnew = this.feSpace.secondMember(l);
                    case 'double'
                        lnew = LinearForm(this.feSpace,l);
                end
                C = plus@AbstractMatrix(this,lnew);
            elseif isequal(class(l),'LinearForm')
                C = plus(l,this);
            end
        end
        function[C] = mtimes(this,u)
            if isa(u,'FE_func')
                C = this.eval(u);
            elseif and(isa(u,'double'),isscalar(u))
                C = mtimes@AbstractMatrix(u,this);
            elseif isa(u,'LinearForm')
                C = mtimes(u,this);
            else
                % The result is no longer a linear form, but still an
                % abstract matrix.
                assert(size(this,2)==size(u,1));
                C = AbstractMatrix(this)*AbstractMatrix(u);
            end
        end
    end
end

















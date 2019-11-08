classdef func
    % Function defined on a mesh, and possibly formed from two
    % components, a FE part which is a FE_func object
    % (then the mesh is that of the FE part), and a
    % handle part, which is simply a function handle defined on the plane.
    
    properties
        feSpace;
        fePart@double;
        handlePart;
    end
    methods
        function[this] = func(arg1,arg2,arg3)
            if nargin == 0
            elseif nargin==1
                switch(class(arg1))
                    case 'FEspace'
                        this = func(arg1,sparse(arg1.ndof,1),[]);
                    case 'func'
                        this = arg1;
                    case 'FE_func',
                        this = func(arg1.feSpace,arg1.fePart,arg1.handlePart);
                    case 'R2toRfunc'
                        this = func(arg1.feSpace,arg1.fePart,arg1.handlePart);
                    case 'double'
                        this = func(R2toRfunc(arg1));
                end
            elseif nargin==2
                if and(isa(arg1,'FE_func'),isa(arg2,'R2toRfunc'))
                    this.feSpace = arg1.feSpace;
                    this.fePart = arg1.fePart;
                    this.handlePart = arg2.handlePart;
                elseif and(isa(arg1,'FEspace'),isa(arg2,'func'))
                    this.feSpace = arg1;
                    if isempty(arg2.fePart)
                        this.fePart = sparse(arg1.ndof,1);
                    else
                        this.fePart = arg2.fePart;
                        assert(length(arg2.fePart)==arg1.ndof);
                    end
                    this.handlePart = arg2.handlePart;
                else
                    error('unexpected...');
                end
            else
                if ~isempty(arg1)
                    if isempty(arg2)
                        arg2 = sparse(arg1.ndof,1);
                    end
                end
                if isempty(arg3)
                    arg3 = @(Z)(0*Z(:,1));
                end
                this.feSpace = arg1;
                this.fePart = arg2;
                this.handlePart = arg3;
            end
        end
        function[] = plot(this,varargin)
            Vh = this.feSpace;
            mesh = Vh.mesh;
            phi = this.feSpace.phi;
            sVec = this.feSpace.sVec;
            M = mesh.M(sVec);
            if ~isempty(this.fePart)
                val1 = phi*this.fePart;
            else
                val1 = 0;
            end
            if ~isempty(this.handlePart)
                val2 = this.handlePart(M);
            else
                val2 = 0;
            end
            plot(sVec,real(val1 + val2),varargin{:});
            if max(abs(imag(val1+val2)))>1e-7    
                warning('Imaginary part of the data ignored');
            end
        end
        function[out] = vals(this)
            Vh = this.feSpace;
            mesh = Vh.mesh;
            phi = this.feSpace.phi;
            sVec = this.feSpace.sVec;
            M = mesh.M(sVec);
            if ~isempty(this.fePart)
                val1 = phi*this.fePart;
            else
                val1 = 0;
            end
            if ~isempty(this.handlePart)
                val2 = this.handlePart(M);
            else
                val2 = 0;
            end
            out = val1 + val2;
        end
        function[uh,f] = parts(this)
            if isempty(this.feSpace)
                uh = [];
            else
                uh = FE_func(this.feSpace,this.fePart);
            end
            
            f = R2toRfunc(this.handlePart);
        end
        function[Vh] = commonSpace(A,B)
            if ~isprop(A,'feSpace')
                if isprop(B,'feSpace')
                    Vh = B.feSpace;
                elseif ~isprop(B,'feSpace')
                    Vh = [];
                end
            elseif ~isprop(B,'feSpace')
                Vh = commonSpace(B,A);
            else
                if and(isempty(A.feSpace),isempty(B.feSpace))
                    Vh = [];
                elseif isempty(A.feSpace)
                    Vh = B.feSpace;
                elseif isempty(B.feSpace)
                    Vh = A.feSpace;
                else
                    Vh = A.feSpace;
                    assert(isequal(B.feSpace,Vh),'arguments do not belong to the same FE space');
                end
                
            end
            
        end
        function[this] = set.feSpace(this,Vh)
            assert(or(isa(Vh,'FEspace'),isempty(Vh)));
            this.feSpace = Vh;
        end
        function[this] = set.handlePart(this,f)
            if isa(this,'FE_func')
                error('The handle part of a ''FE_func'' object must stay unchanged equal to R2toRfunc(0)')
            end
            this.handlePart = f;
        end
        function[A_new,B_new] = operationInterpreter(A,B)
            if isequal(class(A),'func')
                A_new = A;
                switch class(B)
                    case 'func'
                        B_new = B;
                        Vh = commonSpace(A_new,B_new);
                        A_new.feSpace = Vh;
                        B_new.feSpace = Vh;
                    case 'FE_func'
                        B_new = func(B);
                        Vh = commonSpace(A_new,B_new);
                        A_new.feSpace = Vh;
                    case 'R2toRfunc'
                        B_new = func(B);
                        B_new.feSpace = A_new.feSpace;
                    case 'double'
                        B_new = func(B);
                        B_new.feSpace = A_new.feSpace;
                    otherwise
                        error('Incompatible classes %s and %s',...
                            class(A),class(B));
                end
            elseif isequal(class(A),'FE_func')
                switch class(B)
                    case 'FE_func'
                        A_new = A;
                        B_new = B;
                        B_new.feSpace = commonSpace(A,B);
                    case 'R2toRfunc'
                        A_new = func(A);
                        B_new = func(B);
                        B_new.feSpace = A.feSpace;
                        B_new.fePart = sparse(A.feSpace.ndof,1);
                    case 'double'
                        A_new = A;
                        B_new = FE_func(A.feSpace,B);
                    otherwise
                        [B_new,A_new] = operationInterpreter(B,A);
                end
            elseif isequal(class(A),'R2toRfunc')
                switch class(B)
                    case 'R2toRfunc'
                        A_new = A;
                        B_new = B;
                    case 'double'
                        A_new = A;
                        B_new = R2toRfunc(B);
                    otherwise
                        [B_new,A_new] = operationInterpreter(B,A);
                end
            elseif isequal(class(A),'double')
                [B_new,A_new] = operationInterpreter(B,A);
            end
        end
        function[C] = plus(A,B)
            [A_new,B_new] = operationInterpreter(A,B);
            cl = class(A_new);
            assert(isequal(cl,'func'),'unexpected case');
            [uh1,f1] = parts(A_new);
            [uh2,f2] = parts(B_new);
            uh = uh1 + uh2;
            f = f1 + f2;
            C = func(uh,f);
        end
        function[C] = mtimes(A,B)
            if isa(A,'func')
                assert(and(isa(B,'double'),isscalar(B)));
                [uh1,f1] = parts(A);
                uh = B*uh1;
                f = B*f1;
                C = func(uh,f);
            else
                C = mtimes(B,A);
            end
        end
        function[C] = conj(A)
            [uh1,f1] = parts(A);
            uh = conj(uh1);
            f = conj(f1);
            C = func(uh,f);
        end
        function[C] = uminus(A)
            C = -1*A;
        end
        function[C] = minus(A,B)
            C = A + -B;
        end
        function[r] = scal(A,B,Vh)
            
            [A_new,B_new] = operationInterpreter(A,B);
            if or(nargin==2,isempty(Vh))
                Vh = commonSpace(A,B);
                assert(~isempty(Vh),...
                    'You must specify on which Fe space this scalar product must be computed...');
            end
            q = Vh.quad.num;
            [uh1,f1] = parts(A_new);
            [uh2,f2] = parts(B_new);
            M = Vh.Mass;
            l1 = Vh.secondMember(f1,'quadrature',q);
            l2 = Vh.secondMember(f2,'quadrature',q);
            r = (M*uh1|uh2) + (l1|uh2) + (l2|uh1) + Vh.integral(f1*f2);
        end
        function[r] = or(A,B)
            r = scal(A,B,[]);
        end
        
        
        function[] = showLog(this)
            Vh = this.feSpace;
            mesh = Vh.mesh;
            x = linspace(0,1,100)';
            x(1) = x(1) + 0.0001;
            x(end) = x(end) - 0.0001;
            [U,~,sVec] = Vh.dof2point(x);
            M = mesh.M(sVec);
            semilogy(sVec,abs(U*this.fePart + this.handlePart(M)));
        end
    end
    
end


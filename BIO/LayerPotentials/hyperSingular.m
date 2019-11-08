classdef hyperSingular < BIO
    % weighted hypersingular operator N_omega
    properties
        r=1;
        k;
        ln_reg;
        ln_dx_reg;
    end
    
    methods
        function[this] = hyperSingular(VVh,varargin)
            assert(isa(VVh,'FEspace'));
            p = inputParser;
            p.addOptional('k',0);
            p.addOptional('Xdata',VVh.gaussPoints );
            p.addOptional('Op_opt',{})
            p.addOptional('r',1);
            p.parse(varargin{:});
            k = p.Results.k; XX = p.Results.Xdata; AopOpt = p.Results.Op_opt;
            r = p.Results.r;
            this.r = r;
            this.k = k;
            this.Vh = VVh;
            this.ln_reg = this.Vh.regularize(XX,'ln','correcMethod','constantTerm');
            [Mx,My] = this.Vh.regularize(XX,'ln_dx','correcMethod','constantTerm');
            this.ln_dx_reg{1} = Mx;
            this.ln_dx_reg{2} = My;
            if k > 0
                kern = 1i/4*H0Kernel(k);
            else
                kern = (-1/(2*pi))*LogKernel(r);
            end
            this.kernel = kern;
            this.Aop = Op(XX,kern,VVh.gaussPoints,AopOpt{:});
            this.AopOpt = AopOpt;
            this.X = XX;
            this.V = 'V';
        end
        function[] = set_X(this,XX,varargin)
            %             if isequal(XX,this.X)
            %                 return
            %             end
            this.Aop = this.Aop.update_X(XX,varargin{:});
            this.ln_omega_reg = this.Vh.regularizeHypersingular(XX);
            this.X = XX;
        end
        function[this] = remesh(this,N)
            VVh = this.Vh.remesh(N);
            if isequal(this.X,this.Vh.gaussPoints)
                XX = VVh.gaussPoints;
            else
                XX = this.X;
            end
            this = singleLayer(this.k,VVh,XX,this.AopOpt,this.r);
        end
        function[M] = Mat(this)
            % This retunrs the matrix such that (M*U)_i = Hu(x_i)
            Mkern = this.Aop;
            
            M1 = Mkern*(AbstractMatrix.spdiag(this.Vh.W)*this.Vh.dphi) + ...
                -1/(2*pi)*this.ln_dx_reg;
            M2 = Mkern*(AbstractMatrix.spdiag(this.Vh.W)*this.Vh.phi) + ...
                -1/(2*pi)*this.ln_reg;
            M = M1 - this.k^2*M2;
        end
        function[bili] = galerkine(this)
            N = this.Vh.normVecGauss;
            Nx = N(:,1);
            Ny = N(:,2);
            Mkern = this.Aop;
            U = this.Vh.phi; dU = this.Vh.dphi; 
            W = AbstractMatrix.spdiag(this.Vh.W);
            Wx = AbstractMatrix.spdiag(this.Vh.W.*Nx);
            Wy = AbstractMatrix.spdiag(this.Vh.W.*Ny);
            M1 = dU'*W*(Mkern*(W*dU) + ...
                -1/(2*pi)*this.ln_dx_reg);
            M2 = U'*Wx*(Mkern*(Wx*U) + ...
                -1/(2*pi)*this.ln_reg);
            bili = BilinearForm(this.Vh,this.Vh,M1 - this.k^2*M2);
        end
    end
end



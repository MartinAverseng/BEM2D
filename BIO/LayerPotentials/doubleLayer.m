classdef doubleLayer < BIO
    
    properties
        r=1;
        k;
        nscalrRm1_reg;
        correctionMethod = 'constantTerm';
    end
    
    methods
        function[this] = doubleLayer(VVh,varargin)
            p = inputParser;
            p.addOptional('k',0);
            p.addOptional('Xdata',VVh.gaussPoints );
            p.addOptional('Op_opt',{})
            p.addOptional('r',1);
            p.addOptional('correcMethod','constantTerm');
            p.parse(varargin{:});
            k = p.Results.k; XX = p.Results.Xdata; AopOpt = p.Results.Op_opt;
            r = p.Results.r; correcMethod = p.Results.correcMethod;
            this.r = r;
            this.Vh = VVh;
            this.nscalrRm1_reg = this.Vh.regularize(XX,'nscalrRm1','correcMethod',correcMethod);
            if k > 0
                kern = 1i/4*H0Kernel(k);
            else
                kern = (-1/(2*pi))*LogKernel(r);
            end
            this.k = k;
            this.kernel = kern;
            this.Aop = dOp(XX,kern,VVh.gaussPoints,AopOpt{:});
            this.AopOpt = AopOpt;
            this.X = XX;
            this.V = 'V';
        end
        function[this] = set_X(this,XX,varargin)
%             if isequal(XX,this.X)
%                 return
%             end
              this = doubleLayer(this.Vh,'Xdata',XX,varargin{:});
        end
        function[this] = remesh(this,N)
            VVh = this.Vh.remesh(N);
            if isequal(this.X,this.Vh.gaussPoints)
                XX = VVh.gaussPoints;
            else
                XX = this.X;
            end
            this = doubleLayer(this.k,VVh,XX,this.AopOpt,this.r);
        end
        function[M] = Mat(this)
            % This retunrs the matrix such that (M*U)_i = Su(x_i)
            N = this.Vh.normVecGauss;
            Nx = N(:,1);
            Ny = N(:,2);
            Ax = this.Aop.dx;
            Ay = this.Aop.dy;
            Mx = Ax*(AbstractMatrix.spdiag(this.Vh.W.*Nx))*this.testFunc;
            My = Ay*(AbstractMatrix.spdiag(this.Vh.W.*Ny))*this.testFunc;
            M = Mx + My - 1/(2*pi)*this.nscalrRm1_reg;
        end
    end
end



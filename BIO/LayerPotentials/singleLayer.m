classdef singleLayer < BIO
    
    properties
        r=1;
        k;
        ln_reg;
        correctionMethod = 'constantTerm';
    end
    
    methods
        function[this] = singleLayer(VVh,varargin)
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
            this.ln_reg = this.Vh.regularize(XX,'ln','correcMethod',correcMethod);
            if k > 0
                kern = 1i/4*H0Kernel(k);
            else
                kern = (-1/(2*pi))*LogKernel(r);
            end
            this.k = k;
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
            this.ln_reg = this.Vh.regularize(XX,'ln');
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
            % This retunrs the matrix such that (M*U)_i = Su(x_i)
            Mkern = this.Aop;
            M = Mkern*(AbstractMatrix.spdiag(this.Vh.W)*this.testFunc) + ...
                -1/(2*pi)*this.ln_reg;
        end
    end
end



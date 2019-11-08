classdef BIO < handle
    properties
        Vh@FEspace
        kernel@Kernel
        Aop
        X
        V@char
        AopOpt
    end
    methods
        function[this] = BIO(VVh,kernel,XX,VV,varargin)
            if nargin == 0
                return
            end
            this.Vh = VVh;
            this.kernel = kernel;
            this.X = XX;
            this.V = VV;
            YGauss = VVh.gaussPoints;
            this.Aop = Op(XX,kernel,YGauss,varargin{:});
            this.AopOpt = varargin;
        end
        function[Test] = testFunc(this)
            if isequal(this.V,'V')
                Test = this.Vh.phi;
            elseif isequal(this.V,'dV')
                Test = this.Vh.dphi;
            else
                error('Not recognized string');
            end
        end
        function[M] = Mat(this)
            W = AbstractMatrix.spdiag(this.Vh.W);
            M = this.Aop*(W*this.testFunc);
        end
        function[] = set_X(this,XX,varargin)
            this.Aop = this.Aop.update_X(XX,varargin{:});
            this.X = XX;
        end
        function[this] = remesh(this,N)
            Vhh = this.Vh.remesh(N);
            if isequal(this.X,this.Vh.gaussPoints)
                XX = Vhh.gaussPoints;
            else
                XX = this.X;
            end
            this = BIO(Vhh,this.kernel,XX,this.V,this.AopOpt);
        end
        function[vals] = mtimes(S,lambda)
            assert(S.Vh.contains(lambda));
            vals = S.Mat*lambda.v;
        end
        function[bili] = galerkine(this,Wh,U)
            % Returns the bilinear form 
            % b(uh,vh) = \int_{Wh}\int_{Vh}uh(x)K(x,y)vh(y)dx dy
            if ~isequal(this.X,Wh.gaussPoints)
                this.set_X(Wh.gaussPoints);
            end
            if isequal(U,'U')
                Test_left = Wh.phi;
            elseif isequal(U,'dU')
                Test_left = Wh.dphi;
            else
                error('Not recognized string');
            end
            mat = (Test_left.'*AbstractMatrix.spdiag(Wh.W).')*this.Mat;
            bili = BilinearForm(Wh,this.Vh,mat);
        end
    end
end


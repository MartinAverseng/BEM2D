classdef P1disc < FEcell
    
    properties
    end
    
    methods
        function[this] = P1disc()
            this.order = 1;
            this.Nb = 2;
            this.name = 'Discontinuous piecewise affine (P1) elements';
            this.id = 'P1disc';
        end
        function[phi,dxphi] = basis(~,x)
            assert(isequal(x,x(:)));
            phi = [1-x,x];
            dxphi = [-ones(size(x)),ones(size(x))];
        end
        function [T] = dofIndexes(~,mesh)
            T = [2*(1:mesh.nseg)'-1 2*(1:mesh.nseg)'];
        end
        function [Z] = dofCoords(~,mesh)
            X = mesh.vertices(:,1);
            X = [X X]';
            X = X(:);
            X = [X(2:end);X(1)];
            Y = mesh.vertices(:,2);
            Y = [Y Y]';
            Y = Y(:);
            Y = [Y(2:end);Y(1)];
            Z = [X Y];
        end
        function[f] = singularIntegralDictionnary(this,name)
            switch name
                case 'ln'
                    %
                    f{1} = @(X,A,B)(integralLn(1,X,A,B));
                    f{2} = @(X,A,B)(integralLn(2,X,A,B));
                    
                case 'doubleLayer'
                    f{1} = @(X,A,B,N)(integralDoubleLayer(1,X,A,B,N));
                    f{2} = @(X,A,B,N)(integralDoubleLayer(2,X,A,B,N));
                otherwise
                    error('unable to compute singular integral correction for kernel %s : this kernel is unknown',name);
            end
            function [out] = integralLn(b,X,A,B)
                I0 = I_0(X,A,B);
                I1 = I_1(X,A,B);
                [alpha,beta] = this.phi_b(X,A,B);
                out = alpha{b}.*I1 + beta{b}.*I0;
            end
            function [out] = integralDoubleLayer(b,X,A,B,N)
                J0 = J_0(X,A,B,N);
                J1 = J_1(X,A,B,N);
                [alpha,beta] = this.phi_b(X,A,B);
                out = alpha{b}.*J1 + beta{b}.*J0;
            end
        end
        function[f] = singularIntegralDictionnaryDer(this,name)
            switch name
                case 'ln'
                    f{1} = @(X,A,B)(integralLn(1,X,A,B));
                    f{2} = @(X,A,B)(integralLn(2,X,A,B));
            end
            function [out] = integralLn(b,X,A,B)
                I0 = I_0(X,A,B);
                [alpha,~] = this.phi_b(X,A,B);
                out = alpha{b}.*I0;
            end
        end
        function [alpha,beta] = phi_b(~,X,A,B)
            [~,~,~,~,l,~,a,b,~] = parameters_singInt(X,A,B);
            
            alpha{1} = -1./l;
            beta{1} = b./l;
            
            alpha{2} = 1./l;
            beta{2} = -a./l;
            
        end
        function[cb] = constantTerm(~,b,X,A,B)
            [~,~,~,~,~,~,alpha,beta,~] = parameters_singInt(X,repmat(A,size(X,1),1),repmat(B,size(X,1),1));
            if b==1
                cb = beta./(beta-alpha);
            else
                cb = -alpha./(beta-alpha);
            end
        end
    end
    
end
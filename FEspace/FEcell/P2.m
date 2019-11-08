classdef P2 < FEcell
    
    properties
    end
    
    methods
        function[this] = P2()
            this.order = 2;
            this.Nb = 3;
            this.name = 'Continuous piecewise quadratic (P2) elements';
            this.id = 'P2';
        end
        function[phi,dxphi] = basis(~,x)
            assert(isequal(x,x(:)));
            phi = [2*x.^2 - 3*x + 1,-4*x.^2 + 4*x, 2*(1-x).^2 - 3*(1-x) + 1];
            dxphi = [4*x - 3,-8*x + 4, -4*(1-x) + 3];
        end
        function [T] = dofIndexes(~,mesh)
            ndof = mesh.nseg*2 + (~mesh.curve.closed)*1;
            T = [(1:2:ndof-2)',(2:2:ndof-1)',(3:2:ndof)'];
            if mesh.curve.closed
                T = [T;ndof-1,ndof,1];
            end
            assert(size(T,1)==mesh.nseg);
        end
        function [Z] = dofCoords(this,mesh)
            ndof = mesh.nseg*2 + (~mesh.curve.closed)*1;
            [x1,x2,y1,y2] = mesh.edgesCoords;
            T = this.dofIndexes(mesh);
            Z = zeros(ndof,1);
            Z(T(:,1),1) = x1;
            Z(T(:,1),2) = y1;
            Z(T(:,2),1) = (x1+x2)/2;
            Z(T(:,2),2) = (y1+y2)/2;
            Z(T(:,3),1) = x2;
            Z(T(:,3),2) = y2; % This line is very redundant because it should
            % only change at most 1 value if the curve is not closed (otherwise
            % it is useless).
        end
        function[f] = singularIntegralDictionnary(this,name)
            switch name
                case 'ln'                    
                    f{1} = @(X,A,B)(integralLn(1,X,A,B));
                    f{2} = @(X,A,B)(integralLn(2,X,A,B));
                    f{3} = @(X,A,B)(integralLn(3,X,A,B));                    
                case 'doubleLayer'                    
                    f{1} = @(X,A,B,N)(integralDoubleLayer(1,X,A,B,N));
                    f{2} = @(X,A,B,N)(integralDoubleLayer(2,X,A,B,N));
                    f{3} = @(X,A,B,N)(integralDoubleLayer(3,X,A,B,N));     
                otherwise
                    error('unable to compute singular integral correction for kernel %s : this kernel is unknown',name);
            end
            function[out] = integralLn(b,X,A,B)
                I0 = I_0(X,A,B);
                I1 = I_1(X,A,B);
                I2 = I_2(X,A,B);
                [alpha,beta,gamma] = this.phi_b(X,A,B);
                out = alpha{b}.*I2 + beta{b}.*I1 + gamma{b}.*I0;
            end
            function[out] = integralDoubleLayer(b,X,A,B,N)
                J0 = J_0(X,A,B,N);
                J1 = J_1(X,A,B,N);
                J2 = J_2(X,A,B,N);
                [alpha,beta,gamma] = this.phi_b(X,A,B);
                out = alpha{b}.*J2 + beta{b}.*J1 + gamma{b}.*J0;
            end
            
        end
        function[f] = singularIntegralDictionnaryDer(this,name)
            switch name
                case 'ln'
                    f{1} = @(X,A,B)(integralLn(1,X,A,B));
                    f{2} = @(X,A,B)(integralLn(2,X,A,B));
                    f{3} = @(X,A,B)(integralLn(3,X,A,B));
                case 'doubleLayer'
                    f{1} = @(X,A,B,N)(integralDoubleLayer(1,X,A,B,N));
                    f{2} = @(X,A,B,N)(integralDoubleLayer(2,X,A,B,N));
                    f{3} = @(X,A,B,N)(integralDoubleLayer(3,X,A,B,N));
            end
            function [out] = integralLn(b,X,A,B)
                I0 = I_0(X,A,B); 
                I1 = I_1(X,A,B);
                [alpha,beta,~] = this.phi_b(X,A,B);
                out = 2*alpha{b}.*I1 + beta{b}.*I0;
            end
            function [out] = integralDoubleLayer(b,X,A,B,N)
                J0 = J_0(X,A,B,N); 
                J1 = J_1(X,A,B,N);
                [alpha,beta,~] = this.phi_b(X,A,B);
                out = 2*alpha{b}.*J1 + beta{b}.*J0;
            end
        end
        function[alpha,beta,gamma] = phi_b(~,X,A,B)
            [~,~,~,~,l,~,a,b,~] = parameters_singInt(X,A,B);
            alpha{1} = 2./l.^2;
            beta{1} = -4*a./l.^2 - 3./l;
            gamma{1} = 1 + 2*a.^2./l.^2 + 3*a./l;
            
            alpha{2} = -4./l.^2;
            beta{2} =  8*a./l.^2 + 4./l;
            gamma{2} = -4*a.^2./l.^2 - 4*a./l;
            
            alpha{3} = 2./l.^2;
            beta{3} = -4*b./l.^2 + 3./l;
            gamma{3} = 1 + 2*b.^2./l.^2 - 3*b./l;
            
        end
    end
end


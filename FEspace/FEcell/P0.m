classdef P0 < FEcell
    
    properties
    end
    
    methods
        function[this] = P0()
            this.order = 0;
            this.Nb = 1;
            this.name = 'Discontinuous piecewise constant (P0) elements';
            this.id = 'P0';
        end
        function[phi,dxphi] = basis(~,x)
            phi = 0*x + 1;
            dxphi = 0*x;
        end
        function [T] = dofIndexes(~,mesh)
            T = (1:mesh.nseg)';
        end
        function [Z] = dofCoords(~,mesh)
            Z = mesh.middle;
        end
        function[cb] = constantTerm(~,~,X,~,~)
            cb = ones(size(X,1),1); % easy ones(easy,1)
        end
        function[f] = singularIntegralDictionnary(~,name)
            switch name
                case 'ln'
                    % We have to return the function (X,A,B) ->
                    % \int{[AB]}\ln|X-Y|dY
                    f{1} = @I_0;        
                case 'doubleLayer'
                    f{1} = @J_0;
                otherwise
                    error('unable to compute singular integral correction for kernel %s : this kernel is unknown',name);
            end
        end
        
        function [I] = singularIntegralDictionnaryDer(~,~)            
            I = error('Not relevant');
        end
    end
end


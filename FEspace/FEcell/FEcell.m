classdef (Abstract) FEcell
    % Object defining a Finite element 'cell', that is the properties of a
    % finite element basis on the reference element.
    % It is an abstract class, because some properties of a FE cell
    % are shared by all the possible choices of specific FE basis.
    properties
        order; % Order of polynomials involved
        Nb; % Number of basis functions
        name; % Name of the FE.
        id;
    end
    methods
        function[] = plot(this)
            plotRef(this);
        end
        function[N] = ndof(this,mesh)
            N = size(this.dofCoords(mesh),1);
        end
        function [U,dsU,L] = dof2points(this,xhat,mesh)
            
            % * U is the matrix that maps the coordinates of a FE function
            % defined on a mesh, to the values of this function on the
            % points in between. More precisely, given the values phi(x1),
            % ... phi(xN) (in a vector P(i) = pĥi(Ai), of the function phi on the N dof A1, ... , AN,
            % the vector v = U*P contains the v(j) = phi(xj) where the
            % points xj are obtained by dividing each segment [Ai Ai+1]
            % into [Ai+xhat(1)*(Ai+1 - Ai),...,Ai + xhat(end)*(Ai+1 - Ai)].
            % The values xhat must lie between 0 and 1, with extremties
            % excluded.
            % * dsU is analog to U but returning the tangential derivative of the
            % function instead of the function itself.
            % * xhat is an array of values in ]0 1[.
            % * mesh is a MeshCurve object
            
            
            assert(isequal(xhat,xhat(:)));
            assert(sum(ismember([0,1],xhat))==0,...
                ['This function does not work yet when '...
                'the extremities belong to the integration points.']);
            % It is harder to deal with the contrary case.
            % Maybe we will fix this sometime.
            n = length(xhat);
            nseg = mesh.nseg;
            nint = n*nseg;
            Ndof = this.ndof(mesh);
            sizeU = [nint,Ndof];
            
            % Pre-allocating for sparse matrix assembling.
            i = zeros(this.Nb*nseg*n,1);
            j = i; val1 = i; val2 = i;
            T = this.dofIndexes(mesh);
            
            % The local basis that must be copied on each segment
            [phi,dxphi] = this.basis(xhat);
            
            % Lengths for the jacobian
            L = mesh.length;
            for k = 1:n
                for b = 1:this.Nb
                    indexes = (k-1)*this.Nb*nseg + (b-1)*nseg + (1:nseg);
                    lines = k:n:nint;
                    cols = T(:,b);
                    values1 = phi(k,b)*ones(nseg,1);
                    values2 = dxphi(k,b)*1./L;
                    i(indexes) = lines;
                    j(indexes) = cols;
                    val1(indexes) = values1;
                    val2(indexes) = values2;
                end
            end
            U = sparse(i,j,val1,sizeU(1),sizeU(2));
            dsU = sparse(i,j,val2,sizeU(1),sizeU(2));
        end
        function[phi_b_x] = basis_phi_b(this,b,x,der)
            assert(isequal(x,x(:)));
            [phi,dphi] = this.basis(x);
            if der
                phi_b_x = dphi(:,b);
            else
                phi_b_x = phi(:,b);
            end
            
        end
        function [] = plotRef(this)
            % plots the reference element with the basis functions
            x = linspace(0,1,50)';
            [phi,dxphi] = this.basis(x);
            
            X = this.localDof;
            plot(X,0*X,'ks','DisplayName','dof','MarkerSize',12);
            hold on
            [~,dXphi] = this.basis(X);
            plot([0,1],[0,0],'k-','HandleVisibility','off');
            plot([0,1],[0,0],'k+','MarkerSize',9,'HandleVisibility','off')
            for b = 1:this.Nb
                if b==1
                    plot(x,phi(:,b),'b.-','DisplayName','Basis functions');
                    plot(x,dxphi(:,b),'r--','DisplayName','Derivatives');
                else
                    plot(x,phi(:,b),'b.-','HandleVisibility','off');
                    plot(x,dxphi(:,b),'r--','HandleVisibility','off');
                end
                
                plot([x(1) x(end)],[phi(1,b) phi(end,b)],'bo','HandleVisibility','off');
                text(X(b),0.95,['\Phi_' sprintf('%d',b) '(x)'],'Color','b');
                
                plot([x(1) x(end)],[dxphi(1,b) dxphi(end,b)],'r*','HandleVisibility','off');
                text(X(b),dXphi(b,b)+0.05,['\Phi_' sprintf('%d',b) '''(x)'],'Color','r');
            end
            legend show
            xlim([-0.5, 1.5]);
            yl = ylim;
            Delta = diff(yl);
            yl = yl + 0.1*Delta*[-1,2];
            ylim(yl)
            title(this.name);
            set(gca,'XTick',unique([0;X;1]));
        end
        function[X] = localDof(this)
            % Returns the location of the degrees of freedom on the
            % reference element.
            
            % We create a fake curve.
            curve = openline(0,1);
            % We create a mesh with only one segment
            mesh = MeshCurve(curve,1);
            % We retrieve the coordinates of the dofs defined for this
            % specific FE
            Z = this.dofCoords(mesh);
            % Only the x component matters.
            X = Z(:,1);
        end
    end
    methods (Abstract)
        [phi,dxphi] = basis(this,x)
        % Basis functions of the Finite Element and their derivative
        
        
        
        T = dofIndexes(this,mesh)
        % T is a table of size [Nsegment, Nb], and T(i,:) contains the list
        % of the indexes of the dofs contained in this segment.
        
        Z = dofCoords(this,mesh)
        % Returns the coordinates of the dofs for the mesh given in argumen
        
        I = singularIntegralDictionnary(this,name)
        % Singular integrals needed to assemble Bilinear forms.
        
        cb = constantTerm(this,b,X,A,B)
        % cb = constantTerm(b,X,A,B)
        % -X, A and B are three lists of M points in R^2.
        % b is an integer between 1 and Nb.
        % - the output cb is the vector of size M such that, for each
        % k in {1,..,M}, cb(k) is the constant such that if phi(t) is the
        % b-th basis on [A,B] and if t0 = u.XkA, (u = unit vector of AB)
        % then phi - c(Xk) is the restriction of a polynomial vanishing at
        % t=t0.
        
        
        I = singularIntegralDictionnaryDer(this,name)
        
        
        
    end
end


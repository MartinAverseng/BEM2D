classdef FEspace < handle
    % A FE space is an object defining an ensemble of functions defined on
    % a mesh and of the form f(x) = \sum \alpha_i \phi_i(x), where the
    % phi_i are basis functions described by the fe 'cell'.
    
    
    properties
        mesh@MeshCurve;
        fe_cell@FEcell;
        quad@NumericalQuadrature
        phi
        dphi
        sVec
        gaussPoints
        W
        Mass@BilinearForm
        dMass@BilinearForm
        selectedTest
    end
    
    methods
        % Constructor
        function[this] = FEspace(m,c,quad)
            assert(ischar(c));
            switch c
                case 'P0'
                    this.fe_cell = P0;
                case 'P1'
                    this.fe_cell = P1;
                case 'P1disc'
                    this.fe_cell = P1disc;
                case 'P2'
                    this.fe_cell = P2;
                otherwise
                    error('unknown type of finite element space %s, or not supported yet',c);
            end
            this.mesh = m;
            if ~exist('quad','var')||isempty(quad)
                this.quad = inputOrDefaultQuad();
            else
                this.quad = inputOrDefaultQuad(quad);
            end
            
            [U,dU,ssVec,gaussPoints] = dof2point(this.mesh,this.fe_cell,this.quad.xhat);
            this.phi = U;
            this.dphi = dU;
            this.sVec = ssVec;
            this.gaussPoints = gaussPoints;
            this.W = this.quad.weights(this.mesh);
            Wdiag = AbstractMatrix.spdiag(this.W);
            M_mat = AbstractMatrix(U'*Wdiag*U);
            dM_mat = AbstractMatrix(dU'*Wdiag*dU);
            this.Mass = BilinearForm(this,this,M_mat);
            this.dMass = BilinearForm(this,this,dM_mat);
            this.selectedTest = this.phi;
        end
        function[this] = remesh(this,N)
            m = this.mesh.remesh(N);
            this = FEspace(m,this.fe_cell.id,this.quad.num);
            
        end
        function[bool] = contains(this,uh)
            if isa(uh,'FE_func')
                bool = isequal(uh.feSpace,this);
            else
                bool = false;
            end
        end
        
        function[N] = ndof(this)
            % number of degrees of freedom
            N = size(this.phi,2);
        end
        function[N] = nint(this)
            N = size(this.phi,1);
        end
        function[N] = nseg(this)
            % number of segments in the mesh
            N = this.mesh.nseg;
        end
        function[N] = nvert(this)
            % number of vertices in the mesh
            N = this.mesh.nvert;
        end
        function[T] = dofIndexes(this)
            T = this.fe_cell.dofIndexes(this.mesh);
        end
        function[Z] = dofCoords(this)
            Z = this.fe_cell.dofCoords(this.mesh);
        end
        function[N] = normVecGauss(this)
            N = this.mesh.normalVector;
            Nx = N(:,1); Ny = N(:,2);
            Nx = repmat(Nx,1,this.quad.num)'; Ny = repmat(Ny,1,this.quad.num)';
            Nx = Nx(:); Ny = Ny(:); 
            N_norm = sqrt(Nx.^2 + Ny.^2);
            N = [Nx./N_norm, Ny./N_norm];
            
        end
        function[fe_func] = Pi_h(this,func)
            % Creates the vector ndof x 1 of the coordinates of func on the
            % finite element space.
            % func must accept N x 2 arguments and return a N x 1 array
            Z = this.dofCoords;
            assert(isequal(size(Z),[this.ndof,2]));
            assert(isa(func,'R2toRfunc'),sprintf(...
                'Second argument should be of type ''R2toRfunc'' but received %s instead '...
                ,class(func)));
            v = func(Z);
            fe_func = FE_func(this,v);
        end
        function[out] = integral(this,u)
            u = func(u);
            if isempty(u.feSpace)
                u.feSpace = this;
            else
                assert(isequal(this,u.feSpace));
            end
            Wdiag = AbstractMatrix.spdiag(this.W);
            Phi = this.phi;
            gP = this.gaussPoints;
            [uh,f] = parts(u);
            out = sum(Wdiag*Phi*uh.v + Wdiag*f(gP));
        end
        function[l] = secondMember(this,k_func,varargin)
            k_func = func(this,k_func);
            [uh,f] = parts(k_func);
            vals = (f(this.gaussPoints).'.*this.W.')*this.phi;
            M = this.Mass;
            l = LinearForm(this,vals) + M*uh;
        end
        function[l] = Pi_h_normalDerivative(this,fx,fy)
            fx = func(this,fx);
            [uhx,fx] = parts(fx);
            fy = func(this,fy);
            [uhy,fy] = parts(fy);
            assert(sum(abs(uhx.v))==0); 
            assert(sum(abs(uhy.v))==0);
            N = this.mesh.normalVector;
            Nx = N(:,1); Ny = N(:,2);
            valsx = (fx(this.dofCoords).'.*Nx.');
            valsy = (fy(this.dofCoords).'.*Ny.');
            l = FE_func(this,valsx' + valsy');
        end
        function[l] = normalDerivative(this,fx,fy)
            fx = func(this,fx);
            [uhx,fx] = parts(fx);
            fy = func(this,fy);
            [uhy,fy] = parts(fy);
            assert(sum(abs(uhx.v))==0); 
            assert(sum(abs(uhy.v))==0);
            N = this.normVecGauss;
            Nx = N(:,1); Ny = N(:,2);
            valsx = (fx(this.gaussPoints).'.*Nx.'.*this.W.')*this.phi;
            valsy = (fy(this.gaussPoints).'.*Ny.'.*this.W.')*this.phi;
            l = LinearForm(this,valsx + valsy);
        end
        function[a] = meansBasis(this)
            a = sum(spdiags(this.W,0,length(this.W),length(this.W))*this.phi,1)';
        end
        function[dofs] = dofsOnSegments(this,segs)
            T = this.dofIndexes;
            dofs = T(segs,:);
            dofs = unique(dofs(:));
        end
        function[I] = dofCloseToPoints(this,X,threshold)
            I = isClose(X,this.dofCoords,threshold);
        end
        function[Mat] = regularize(this,X,singKernel_id,varargin)
            % Creates the Matrix  Mat such that Mat*U = \int_{\Gamma} kernel(\abs{X-y})u(y)dy -
            % \intApprox{\Gamma} kernel(\abs{X-y})u(y)dy
            % where \intApprox would be what one obtains with a simple numerical
            % quadrature, and u = sum_{i}U_i phi_i(y)
            % (this matrix is the corrective part accounting for singularities).
            p = inputParser;
            p.addOptional('derivative',false);
            p.addOptional('correcMethod','constantTerm');
            p.parse(varargin{:});
            if strcmp(p.Results.correcMethod,'constantTerm')
                Mat = this.correctionConstantTerm(X,singKernel_id);
            end
        end
        function[Mat] = correctionConstantTerm(this,X,singKernel_id)
            % Vh a weighted FE space.
            % This function computes the sparse matrix
            % Aik = Cik(\int_[this]K(X_k,y)  - \intApprox_[this]K(X_k,y))
            % where the non-zeros indexes are chosen when X_k is close to
            % the segments attached to the dof number i. Cik is the
            % piecewise constant function on each segment of the mesh such
            % that the i-th basis function can be written on each segment
            % as phi_i(Y) = Cik(Y) + someConstant times (X_k - Y)
            
            switch singKernel_id
                case 'ln'
                    Mat = correcLn(this,X);
                case 'ln_dx'
                    Mat = correcLnDx(this,X);
                case 'nscalrRm1'
                    Mat = correcrRm1(this,X);
            end
        end
        function [Mat] =  correcLn(this,X)
            
            
            M = size(X,1); % Number of points in X
            me = this.mesh;
            l = me.length; % vector containing length of each segment in the mesh.
            ndof = this.ndof;
            T = this.dofIndexes;
            feCell = this.fe_cell;
            Nb = feCell.Nb;
            singK = logSingK;
            func = singK.k.func;
            
            EDG = this.mesh.edgesCoords; % edges [A, B] with A = [Ax, Ay], B = [Bx, By]
            Ax = EDG(:,1); Ay = EDG(:,2); Bx = EDG(:,3); By = EDG(:,4);
            lines = cell(Nb,1); cols = cell(Nb,1); vals = cell(Nb,1);
            threshold = 1.1*max(l);
            I = isClose(X,this.dofCoords,threshold);
            dofsToLoop = find(~cellfun('isempty',I))';
            for i = dofsToLoop %dofsToLoop
                for b = 1:Nb
                    segNum = find(T(:,b)==i); % On which segment is
                    if ~isempty(segNum)
                        
                        lines_k = I{i}';
                        X_k = X(lines_k,:);
                        % the b-th basis function attached to dof i.
                        %                     X_k = X(I{i}',:);
                        A = [Ax(segNum) Ay(segNum)]; B = [Bx(segNum) By(segNum)]; % the segment is [A,B].
                        
                        
                        Mid = (A + B)/2; % segment middle.
                        lAB2 = l(segNum)^2; % squared length of the segment.
                        
                        
                        
                        % We select in X the points that are within three balls centered
                        % respectively on A,B and Mid, and of radius lAB.
                        AX2 = (X_k(:,1) - Ax(segNum)).^2 + (X_k(:,2) - Ay(segNum)).^2;
                        BX2 = (X_k(:,1) - Bx(segNum)).^2 + (X_k(:,2) - By(segNum)).^2;
                        MX2 = (X_k(:,1) - Mid(1)).^2 + (X_k(:,2) - Mid(2)).^2;
                        ks = find(or(AX2 < lAB2, or(BX2 < lAB2,MX2<lAB2)));
                        lines_k = lines_k(ks);
                        X_k = X_k(ks,:); % seclected points in X.
                        %for each k, parameter cb(k) such that phi_b(Y) = cb(k) + C*(X_k(k,:) - Y)
                        % where C is some constant.
                        approxInt = this.I0approx(func,X_k,A,B,segNum);
                        cb = this.fe_cell.constantTerm(b,X_k,A,B);
                        
                        % Compute the approximation that was used for
                        % \int_{[A,B]}G(X_k(k,:),Y) dY for each k.
                        
                        
                        % Compute the exact integral to replace :
                        exactInt = this.I0exact(singK, X_k,A,B,segNum); % method that computes
                        % exactly \int_{[A,B]} G(X_k(k,:),Y) dY for each k.
                        
                        % Store the correction.
                        correc_kib = cb.*(exactInt - approxInt);
                        lines{b} = [lines{b};lines_k];
                        cols{b} = [cols{b};0*ks + i];
                        vals{b} = [vals{b};correc_kib];
                    end
                end
            end
            % We have a corrective patch for each b on the whole mesh. Sum
            % them to get the final correction.
            Mat = sparse(lines{1},cols{1},vals{1},M,ndof);
            for b = 2:Nb
                Mat = Mat + sparse(lines{b},cols{b},vals{b},M,ndof);
            end
        end
        function [Matx,Maty] =  correcLnDx(this,X)
            
            N = this.mesh.normalVector;
            M = size(X,1); % Number of points in X
            me = this.mesh;
            l = me.length; % vector containing length of each segment in the mesh.
            ndof = this.ndof;
            T = this.dofIndexes;
            feCell = this.fe_cell;
            Nb = feCell.Nb;
            singK = logSingK;
            func = singK.k.func;
            
            EDG = this.mesh.edgesCoords; % edges [A, B] with A = [Ax, Ay], B = [Bx, By]
            Ax = EDG(:,1); Ay = EDG(:,2); Bx = EDG(:,3); By = EDG(:,4);
            lines = cell(Nb,1); cols = cell(Nb,1); vals1 = cell(Nb,1); vals2 = cell(Nb,1);
            threshold = 1.1*max(l);
            I = isClose(X,this.dofCoords,threshold);
            dofsToLoop = find(~cellfun('isempty',I))';
            for i = dofsToLoop %dofsToLoop
                for b = 1:Nb
                    segNum = find(T(:,b)==i); % On which segment is
                    if ~isempty(segNum)
                        
                        lines_k = I{i}';
                        X_k = X(lines_k,:);
                        % the b-th basis function attached to dof i.
                        %                     X_k = X(I{i}',:);
                        A = [Ax(segNum) Ay(segNum)]; B = [Bx(segNum) By(segNum)]; % the segment is [A,B].
                        Nx_i = N(segNum,1);
                        Ny_i = N(segNum,1);
                        Mid = (A + B)/2; % segment middle.
                        lAB2 = l(segNum)^2; % squared length of the segment.
                        
                        
                        
                        % We select in X the points that are within three balls centered
                        % respectively on A,B and Mid, and of radius lAB.
                        AX2 = (X_k(:,1) - Ax(segNum)).^2 + (X_k(:,2) - Ay(segNum)).^2;
                        BX2 = (X_k(:,1) - Bx(segNum)).^2 + (X_k(:,2) - By(segNum)).^2;
                        MX2 = (X_k(:,1) - Mid(1)).^2 + (X_k(:,2) - Mid(2)).^2;
                        ks = find(or(AX2 < lAB2, or(BX2 < lAB2,MX2<lAB2)));
                        lines_k = lines_k(ks);
                        X_k = X_k(ks,:); % seclected points in X.
                        %for each k, parameter cb(k) such that phi_b(Y) = cb(k) + C*(X_k(k,:) - Y)
                        % where C is some constant.
                        approxInt = this.I0approx(func,X_k,A,B,segNum);
                        cb = this.fe_cell.linearTerm(b,X_k,A,B);
                        
                        % Compute the approximation that was used for
                        % \int_{[A,B]}G(X_k(k,:),Y) dY for each k.
                        
                        
                        % Compute the exact integral to replace :
                        exactInt = this.I0exact(singK, X_k,A,B,segNum); % method that computes
                        % exactly \int_{[A,B]} G(X_k(k,:),Y) dY for each k.
                        
                        % Store the correction.
                        correc_kib = cb.*(exactInt - approxInt);
                        lines{b} = [lines{b};lines_k];
                        cols{b} = [cols{b};0*ks + i];
                        vals1{b} = [vals{b};correc_kib*Nx_i];
                        vals2{b} = [vals{b};correc_kib*Ny_i];
                    end
                end
            end
            % We have a corrective patch for each b on the whole mesh. Sum
            % them to get the final correction.
            Matx = sparse(lines{1},cols{1},vals1{1},M,ndof);
            Maty = sparse(lines{1},cols{1},vals2{1},M,ndof);
            for b = 2:Nb
                Matx = Matx + sparse(lines{b},cols{b},vals{b},M,ndof);
                Maty = Maty + sparse(lines{b},cols{b},vals{b},M,ndof);
            end
        end
        function[Mat]= correcrRm1(this,X)
            M = size(X,1); % Number of points in X
            me = this.mesh;
            l = me.length; % vector containing length of each segment in the mesh.
            ndof = this.ndof;
            normalV = this.mesh.normalVector;
            normalV = normalV./sqrt(normalV(:,1).^2 + normalV(:,2).^2); 
            T = this.dofIndexes;
            feCell = this.fe_cell;
            Nb = feCell.Nb;
            
            EDG = this.mesh.edgesCoords; % edges [A, B] with A = [Ax, Ay], B = [Bx, By]
            Ax = EDG(:,1); Ay = EDG(:,2); Bx = EDG(:,3); By = EDG(:,4);
            lines = cell(Nb,1); cols = cell(Nb,1); vals = cell(Nb,1);
            threshold = 1.5*max(l);
            I = isClose(X,this.dofCoords,threshold);
            dofsToLoop = find(~cellfun('isempty',I))';
            for i = dofsToLoop %dofsToLoop
                for b = 1:Nb
                    segNum = find(T(:,b)==i); % On which segment is
                    if ~isempty(segNum)
                        nn = normalV(segNum,:);
                        lines_k = I{i}';
                        X_k = X(lines_k,:);
                        % the b-th basis function attached to dof i.
                        %                     X_k = X(I{i}',:);
                        A = [Ax(segNum) Ay(segNum)]; B = [Bx(segNum) By(segNum)]; % the segment is [A,B].
                        
                        
                        Mid = (A + B)/2; % segment middle.
                        lAB2 = l(segNum)^2; % squared length of the segment.
                        
                        
                        
                        % We select in X the points that are within three balls centered
                        % respectively on A,B and Mid, and of radius lAB.
                        AX2 = (X_k(:,1) - Ax(segNum)).^2 + (X_k(:,2) - Ay(segNum)).^2;
                        BX2 = (X_k(:,1) - Bx(segNum)).^2 + (X_k(:,2) - By(segNum)).^2;
                        MX2 = (X_k(:,1) - Mid(1)).^2 + (X_k(:,2) - Mid(2)).^2;
                        ks = find(or(AX2 < lAB2, or(BX2 < lAB2,MX2<lAB2)));
                        lines_k = lines_k(ks);
                        X_k = X_k(ks,:); % seclected points in X.
                        q = this.quad.num;
                        % Gauss points and weights on segment
                        idx = (segNum-1)*q + (1:q)'; %
                        Yq = this.gaussPoints(idx,:); % points
                        Wq = this.W(idx); % weights
                        %for each k, parameter cb(k) such that phi_b(Y) = cb(k) + C*(X_k(k,:) - Y)
                        % where C is some constant.
                        [~,~,~,~,~,~,alpha,beta,~] = parameters_singInt(X_k,A,B);
                        d = (X_k(:,1) - A(1))*nn(1) + (X_k(:,2) - A(2))*nn(2);

                        Nx = size(X_k,1);
                        Ny = size(Yq,1);
                        x1_y1 = repmat(X_k(:,1),1,Ny) - repmat(Yq(:,1)',Nx,1);
                        x2_y2 = repmat(X_k(:,2),1,Ny) - repmat(Yq(:,2)',Nx,1);
                        rxy = sqrt(x1_y1.^2 + x2_y2.^2);
                        rxy(rxy < 1e-15) = 1e-15;
                        M1 = 1./(rxy).*(-x1_y1./rxy);
                        M2 = 1./(rxy).*(-x2_y2./rxy);
                        approxInt = nn(1)*M1*Wq + nn(2)*M2*Wq;
                        
                        cb = this.fe_cell.constantTerm(b,X_k,A,B);
                        
                        % Compute the approximation that was used for
                        % \int_{[A,B]}G(X_k(k,:),Y) dY for each k.
                        
                        
                        % Compute the exact integral to replace :
                        exactInt = (-atan(beta./d) + atan(alpha./d)); 
                        exactInt(abs(d)<1e-12) = 0; % method that computes
                        % exactly \int_{[A,B]} G(X_k(k,:),Y) dY for each k.
                        
                        % Store the correction.
                        correc_kib = cb.*(exactInt - approxInt);
                        lines{b} = [lines{b};lines_k];
                        cols{b} = [cols{b};0*ks + i];
                        vals{b} = [vals{b};correc_kib];
                    end
                end
            end
            % We have a corrective patch for each b on the whole mesh. Sum
            % them to get the final correction.
            Mat = sparse(lines{1},cols{1},vals{1},M,ndof);
            for b = 2:Nb
                Mat = Mat + sparse(lines{b},cols{b},vals{b},M,ndof);
            end
        end
        function[res] = I0approx(this,func,X_k,~,~,segNum)
            q = this.quad.num;
            % Gauss points and weights on segment
            idx = (segNum-1)*q + (1:q)'; %
            Yq = this.gaussPoints(idx,:); % points
            Wq = this.W(idx); % weights
            rXkYq = sqrt((X_k(:,1)-Yq(:,1)').^2+(X_k(:,2) - Yq(:,2)').^2);
            rXkYq(rXkYq < 1e-15) = 1e-15;
            % in the code of Op so we have to correct it coherently
            Gkq = func(rXkYq);
            %             Gkq(or(isnan(Gkq),isinf(Gkq))) = 0; % this is a global convention
            res = Gkq*Wq;
        end
        function[res] = I0exact(~,singK,X_k,A,B,~)
            res = singK.I0seg(X_k,A,B);
        end
        
        
        
        
        
        %% Deprecated functions :
        
        
        
        %         function[I] = segmentsCloseToDofs(this,threshold)
        %             % Maybe to be removed
        %             hmax = max(this.mesh.length);
        %             Middles = this.mesh.middle;
        %             I = isClose(Middles,this.dofCoords,threshold + 1.5*hmax);
        %         end
        %         function[I] = TheirSegmentsCloseToMyDofs(me,them,threshold)
        %             % Maybe to be removed
        %             mymesh = me.mesh;
        %             myDofs = me.dofCoords;
        %             theirmesh = them.mesh;
        %             hmax = max(mymesh.length);
        %             theirMiddles = theirmesh.middles;
        %             I = isClose(theirMiddles,myDofs,threshold + 1.5*hmax);
        %         end
        %         function[I] = segmentsCloseToDofsOf(me,them,threshold)
        %             % Maybe to be removed
        %             mymesh = me.mesh;
        %             myMiddles = mymesh.middle;
        %             theirmesh = them.mesh;
        %             hmax = max(theirmesh.length);
        %             theirDofs = them.dofCoords;
        %             I = isClose(myMiddles,theirDofs,threshold + 1.5*hmax);
        %         end
        
        % Deprecated, to be removed
        %         function[out] = BIO(this,varargin)
        %             %% Input arguments parsing
        %             p = inputParser;
        %             p.addOptional('k','none');
        %             p.addOptional('dk','none');
        %             p.addOptional('varChange',@(Z)(Z)) %
        %             p.addOptional('singularIntegral',NaN);
        %             p.addOptional('singularIntegralDer',NaN);
        %             p.addOptional('thresholdSingInt',NaN);
        %             p.addOptional('refineQuad',20);
        %             p.addOptional('compressSBD',true);
        %             p.parse(varargin{:});
        %             vars = p.Results;
        %             varChange = vars.varChange;
        %             singInt = vars.singularIntegral;
        %             singIntDer = vars.singularIntegralDer;
        %             threshold = vars.thresholdSingInt;
        %             compressSBD = vars.compressSBD;
        %             k = vars.k;
        %             dk = vars.dk;
        %             q = this.quad;
        %             XGauss = this.gaussPoints;
        %             quadRefine = NumericalQuadrature(vars.refineQuad);
        %             numRefine = quadRefine.num;
        %             if isequal(quadRefine.num,q.num)
        %                 Urefine = this.U;
        %                 dUrefine = this.dU;
        %                 Wrefine = this.W;
        %                 Xrefine = XGauss;
        %             else
        %                 [Urefine,dUrefine,~,Xrefine] = ...
        %                     dof2point(this.mesh,this.fe_cell,quadRefine.xhat);
        %                 Wrefine = quadRefine.weights(this.mesh);
        %             end
        %             %% General variables.
        %             Vh = this;
        %             m = Vh.mesh;
        %             refCell = Vh.cell;
        %
        %             WW = this.W;
        %             nint = size(XGauss,1);
        %             nint_refine = size(Xrefine,1);
        %             ndof = Vh.ndof;
        %             T = Vh.dofIndexes;
        %             Nb = refCell.Nb;
        %             [x1,x2,y1,y2] = m.edgesCoords;
        %             UU = this.U;
        %             dUU = this.dU;
        %             num = q.num;
        %             %% Creation of the Kernel matrix
        %
        %             if compressSBD
        %                 tol = 1e-9;
        %                 a = 1/(2*sqrt(nint));
        %                 b = 1;
        %                 if isequal(k,'none')
        %                     Aop = AbstractMatrix(@(x)(0*x),nint,nint);
        %                 else
        %                     Aop = Op(varChange(XGauss),varChange(XGauss),k,a,b,tol);
        %                 end
        %                 if isequal(dk,'none')
        %                     dAop = AbstractMatrix(@(x)(0*x),nint,nint);
        %                 else
        %                     dAop = Op(varChange(XGauss),varChange(XGauss),dk,a,b,tol);
        %                 end
        %             else
        %                 X1 = XGauss(:,1); X2 = XGauss(:,2);
        %                 DIFF1 = repmat(X1,1,nint) - repmat(X1',nint,1);
        %                 DIFF2 = repmat(X2,1,nint) - repmat(X2',nint,1);
        %                 DIFFNorm = sqrt(DIFF1.^2 + DIFF2.^2);
        %                 if isequal(k,'none')
        %                     Aop = AbstractMatrix(@(x)(0*x),nint,nint);
        %                 else
        %                     Mat = k.func(DIFFNorm);
        %                     Mat(or(isnan(Mat),isinf(Mat)))=0;
        %                     Aop = AbstractMatrix(Mat);
        %
        %                 end
        %                 if isequal(dk,'none')
        %                     dAop = AbstractMatrix(@(x)(0*x),nint,nint);
        %                 else
        %                     Mat = dk.func(DIFFNorm);
        %                     Mat(or(isnan(Mat),isinf(Mat)))=0;
        %                     dAop = AbstractMatrix(Mat);
        %                 end
        %             end
        %
        %             %% Correction of singular integrals
        %             correc = sparse([],[],[],nint,ndof);
        %             fineQuad = sparse([],[],[],nint_refine,ndof);
        %
        %             correcDer = sparse([],[],[],nint,ndof);
        %             fineQuadDer = sparse([],[],[],nint_refine,ndof);
        %             % Correction for the non-derivative part (the k part)
        %             if ~isnan2(singInt)
        %                 % * searching where the integrals are singular
        %                 if isnan(threshold)
        %                     % Proximity threshold, under which inegral are considered
        %                     % singular.
        %                     threshold = max(m.length)/2;
        %                 end
        %                 I = Vh.segmentsCloseToDofs(threshold);
        %
        %
        %                 lines = [];
        %                 cols = lines;
        %                 vals = lines;
        %                 lines_refine = [];
        %                 cols_refine = lines_refine;
        %                 vals_refine = lines_refine;
        %                 for dof = 1:ndof
        %                     neighbourSegs = I{dof};
        %                     A1 = repmat((1:num)',1,length(neighbourSegs));
        %                     A2 = repmat((neighbourSegs-1)*num,num,1);
        %                     A1refine = repmat((1:numRefine)',1,length(neighbourSegs));
        %                     A2refine = repmat((neighbourSegs-1)*numRefine,numRefine,1);
        %                     lines_dof = A1(:) + A2(:);
        %                     cols_dof = 0*lines_dof + dof;
        %                     lines_dof_refine = A1refine(:) + A2refine(:);
        %                     cols_dof_refine = 0*lines_dof_refine + dof;
        %                     X_k = XGauss(lines_dof,:);
        %                     X_krefine = Xrefine(lines_dof_refine,:);
        %                     vals_dof_refine = zeros(length(lines_dof_refine),Nb);
        %
        %                     for b = 1:Nb;
        %                         seg = find(T(:,b)==dof);
        %                         assert(length(seg)<=1);
        %                         if isempty(seg)
        %                             vals_dof_refine(:,b) = 0;
        %                         else
        %                             A = repmat([x1(seg) y1(seg)],length(lines_dof_refine),1);
        %                             B = repmat([x2(seg) y2(seg)],length(lines_dof_refine),1);
        %                             vals_dof_refine(:,b) = singInt{b}(X_krefine,A,B);
        %                         end
        %                     end
        %
        %                     % Old values that need to be substracted :
        %                     l = dof;
        %                     phi_l = UU(:,l);
        %                     p_ind = phi_l~=0;
        %                     X_p = XGauss(p_ind,:);
        %                     W_p = WW(p_ind);
        %                     U_p = phi_l(p_ind);
        %                     n_p = nnz(phi_l);
        %                     n_k = length(lines_dof);
        %                     rXkXp = sqrt((repmat(X_k(:,1),1,n_p) - repmat(X_p(:,1)',n_k,1)).^2 ...
        %                         + (repmat(X_k(:,2),1,n_p) - repmat(X_p(:,2)',n_k,1)).^2);
        %                     Gkp = k.func(rXkXp);
        %                     Gkp(or(isinf(Gkp),isnan(Gkp))) = 0;
        %                     valCorrec = Gkp*(W_p.*U_p);
        %
        %                     % Assembling matrix
        %                     lines = [lines;lines_dof]; % Is it possible to know the
        %                     % size of lines, cols and vals in advance ?
        %                     cols = [cols;cols_dof];
        %                     vals = [vals;-valCorrec];
        %
        %                     lines_refine = [lines_refine;lines_dof_refine];
        %                     cols_refine = [cols_refine;cols_dof_refine];
        %                     vals_refine = [vals_refine;sum(vals_dof_refine,2)];
        %
        %                 end
        %
        %                 correc = sparse(lines,cols,vals,nint,ndof);
        %                 fineQuad = sparse(lines_refine,cols_refine,vals_refine,nint_refine,ndof);
        %             end
        %             % Correction for the derivative part (the dk part)
        %             if ~isnan2(singIntDer)
        %                 % * searching where the integrals are singular
        %                 if isnan(threshold)
        %                     % Proximity threshold, under which inegral are considered
        %                     % singular.
        %                     threshold = max(m.length)/2;
        %                 end
        %                 I = Vh.segmentsCloseToDofs(threshold);
        %
        %
        %                 lines = [];
        %                 cols = lines;
        %                 vals = lines;
        %                 lines_refine = [];
        %                 cols_refine = lines_refine;
        %                 vals_refine = lines_refine;
        %                 for dof = 1:ndof
        %                     neighbourSegs = I{dof};
        %                     A1 = repmat((1:num)',1,length(neighbourSegs));
        %                     A2 = repmat((neighbourSegs-1)*num,num,1);
        %                     A1refine = repmat((1:numRefine)',1,length(neighbourSegs));
        %                     A2refine = repmat((neighbourSegs-1)*numRefine,numRefine,1);
        %                     lines_dof = A1(:) + A2(:);
        %                     cols_dof = 0*lines_dof + dof;
        %                     lines_dof_refine = A1refine(:) + A2refine(:);
        %                     cols_dof_refine = 0*lines_dof_refine + dof;
        %                     X_k = XGauss(lines_dof,:);
        %                     X_krefine = Xrefine(lines_dof_refine,:);
        %                     vals_dof_refine = zeros(length(lines_dof_refine),Nb);
        %
        %
        %                     for b = 1:Nb;
        %                         seg = find(T(:,b)==dof);
        %                         assert(length(seg)<=1);
        %                         if isempty(seg)
        %                             vals_dof_refine(:,b) = 0;
        %                         else
        %                             A = repmat([x1(seg) y1(seg)],length(lines_dof_refine),1);
        %                             B = repmat([x2(seg) y2(seg)],length(lines_dof_refine),1);
        %                             vals_dof_refine(:,b) = singIntDer{b}(X_krefine,A,B);
        %                         end
        %                     end
        %
        %                     % Old values that need to be substracted :
        %                     l = dof;
        %                     dphi_l = dUU(:,l);
        %                     p_ind = dphi_l~=0;
        %                     X_p = XGauss(p_ind,:);
        %                     W_p = WW(p_ind);
        %                     dU_p = dphi_l(p_ind);
        %                     n_p = nnz(dphi_l);
        %                     n_k = length(lines_dof);
        %                     rXkXp = sqrt((repmat(X_k(:,1),1,n_p) - repmat(X_p(:,1)',n_k,1)).^2 ...
        %                         + (repmat(X_k(:,2),1,n_p) - repmat(X_p(:,2)',n_k,1)).^2);
        %                     Gkp = dk.func(rXkXp);
        %                     Gkp(rXkXp<1e-10) = dk.func(1e-10);
        %                     valCorrec = Gkp*(W_p.*dU_p);
        %
        %                     % Assembling matrix
        %                     lines = [lines;lines_dof]; % Is it possible to know the
        %                     % size of lines, cols and vals in advance ?
        %                     cols = [cols;cols_dof];
        %                     vals = [vals;-valCorrec];
        %
        %                     lines_refine = [lines_refine;lines_dof_refine];
        %                     cols_refine = [cols_refine;cols_dof_refine];
        %                     vals_refine = [vals_refine;sum(vals_dof_refine,2)];
        %
        %                 end
        %                 correcDer = sparse(lines,cols,vals,nint,ndof);
        %                 fineQuadDer = sparse(lines_refine,cols_refine,vals_refine,nint_refine,ndof);
        %             end
        %             %% Creation of the bilinear form
        %             mv_p = @(x)(...
        %                 real(...
        %                 UU'*(WW.*(Aop.mv_prod(WW.*(UU*x)) + correc*x))...
        %                 +	Urefine'*(Wrefine.*(fineQuad*x)))...
        %                 + real(...
        %                 dUU'*(WW.*(dAop.mv_prod(WW.*(dUU*x)) + correcDer*x))...
        %                 + dUrefine'*(Wrefine.*(fineQuadDer*x)))...
        %                 );
        %             out = BilinearForm(this,mv_p);
        %         end
        % Deprecated, to be removed
        %         function[S] = singleLayerPotential(this,k,varargin)
        %
        %             p = inputParser;
        %             p.addOptional('r',1);
        %             p.addOptional('compressSBD',true);
        %             p.parse(varargin{:});
        %             vars = p.Results;
        %             r = vars.r;
        %             compressSBD = vars.compressSBD;
        %             if nargin == 1
        %                 k=0;
        %             elseif k~=0
        %                 if r~=1
        %                     warning('value of parameter ''r'' ignored because k~=0')
        %                 end
        %                 r = k;
        %             end
        %             assert(r>0);
        %             G0 = LogKernel(1);
        %             singInt = this.fe_cell.singularIntegralDictionnary('ln');
        %             a = this.meansBasis();
        %             S0 = -1/(2*pi)*this.BIO('k',G0,'singularIntegral',singInt,'compressSBD',compressSBD) + -1/(2*pi)*log(r)*AbstractMatrix(@(x)(a*(a'*x)),this.ndof,this.ndof);
        %             if k ~= 0
        %                 Gk_re = HelmholtzPerturb(k);
        %                 Gk_im = J0Kernel(k);
        %                 Sk_re = this.BIO('k',Gk_re,'compressSBD',compressSBD);
        %                 Sk_im = 1/4*this.BIO('k',Gk_im,'compressSBD',compressSBD);
        %                 S = S0 + Sk_re + 1i*Sk_im;
        %             else
        %                 S = S0;
        %             end
        %         end
        % Deprecated, to be removed
        %         function[T] = doubleLayerPotential(this,varargin)
        %             % For now, I don't know how to compress the kernel, so full
        %             % version.
        %             p = inputParser;
        %             p.addOptional('refineQuad',6);
        %             p.addOptional('thresholdSingInt',NaN);
        %             p.parse(varargin{:});
        %             vars = p.Results;
        %             X = this.zVec;
        %             q = this.quad;
        %             threshold = vars.thresholdSingInt;
        %             quadRefine = NumericalQuadrature(vars.refineQuad);
        %             numRefine = quadRefine.num;
        %             if isequal(quadRefine.num,q.num)
        %                 Urefine = this.U;
        %                 Wrefine = this.W;
        %                 Xrefine = X;
        %             else
        %                 [Urefine,~,~,Xrefine] = ...
        %                     dof2point(this.mesh,this.fe_cell,quadRefine.xhat);
        %                 Wrefine = quadRefine.weights(this.mesh);
        %             end
        %             Vh = this;
        %             m = Vh.mesh;
        %             refCell = Vh.cell;
        %             WW = this.W;
        %             nint = size(X,1);
        %             nint_refine = size(Xrefine,1);
        %             ndof = Vh.ndof;
        %             Tab = Vh.dofIndexes;
        %             Nb = refCell.Nb;
        %             [x1,x2,y1,y2] = m.edgesCoords;
        %             UU = this.U;
        %             num = q.num;
        %
        %
        %             X1 = X(:,1);
        %             X2 = X(:,2);
        %             DIFF1 = repmat(X1,1,nint) - repmat(X1',nint,1);
        %             DIFF2 = repmat(X2,1,nint) - repmat(X2',nint,1);
        %             DIFFNorm = sqrt(DIFF1.^2 + DIFF2.^2);
        %             A1 = DIFF1./DIFFNorm.^2;
        %             A2 = DIFF2./DIFFNorm.^2;
        %             A1(or(isnan(A1),isinf(A1))) = 0;
        %             A2(or(isnan(A2),isinf(A2))) = 0;
        %             N = m.normalVector;
        %             N = N./repmat(sqrt(N(:,1).^2 + N(:,2).^2),1,2);
        %             N1 = spdiags(reshape(repmat(N(:,1)',num,1),nint,1),0,nint,nint);
        %             % Replicate the normal vector (first component) for each integration point.
        %             % And make a sparse diag out of it.
        %             N2 = spdiags(reshape(repmat(N(:,2)',num,1),nint,1),0,nint,nint);
        %             Mat = A1*N1+ A2*N2;
        %             Aop = AbstractMatrix(Mat);
        %             %% Correction of singular integrals
        %
        %             singInt = refCell.singularIntegralDictionnary('doubleLayer');
        %
        %             % Correction for the non-derivative part (the k part)
        %
        %             if isnan(threshold)
        %                 % Proximity threshold, under which inegral are considered
        %                 % singular.
        %                 threshold = max(m.length)/2;
        %             end
        %             I = Vh.segmentsCloseToDofs(threshold);
        %
        %             lines = [];
        %             cols = lines;
        %             vals = lines;
        %             lines_refine = [];
        %             cols_refine = lines_refine;
        %             vals_refine = lines_refine;
        %             for dof = 1:ndof
        %                 neighbourSegs = I{dof};
        %                 A1 = repmat((1:num)',1,length(neighbourSegs));
        %                 A2 = repmat((neighbourSegs-1)*num,num,1);
        %                 A1refine = repmat((1:numRefine)',1,length(neighbourSegs));
        %                 A2refine = repmat((neighbourSegs-1)*numRefine,numRefine,1);
        %                 lines_dof = A1(:) + A2(:);
        %                 cols_dof = 0*lines_dof + dof;
        %                 lines_dof_refine = A1refine(:) + A2refine(:);
        %                 cols_dof_refine = 0*lines_dof_refine + dof;
        %                 X_krefine = Xrefine(lines_dof_refine,:);
        %                 vals_dof_refine = zeros(length(lines_dof_refine),Nb);
        %
        %                 for b = 1:Nb;
        %                     seg = find(Tab(:,b)==dof);
        %                     assert(length(seg)<=1);
        %                     if isempty(seg)
        %                         vals_dof_refine(:,b) = 0;
        %                     else
        %                         A = repmat([x1(seg) y1(seg)],length(lines_dof_refine),1);
        %                         B = repmat([x2(seg) y2(seg)],length(lines_dof_refine),1);
        %                         Nloc = repmat(N(seg,:),length(lines_dof_refine),1);
        %                         vals_dof_refine(:,b) = singInt{b}(X_krefine,A, B,Nloc);
        %                     end
        %                 end
        %
        %                 % Old va%         function[T] = doubleLayerPotential(this,varargin)
        %             % For now, I don't know how to compress the kernel, so full
        %             % version.
        %             p = inputParser;
        %             p.addOptional('refineQuad',6);
        %             p.addOptional('thresholdSingInt',NaN);
        %             p.parse(varargin{:});
        %             vars = p.Results;
        %             X = this.zVec;
        %             q = this.quad;
        %             threshold = vars.thresholdSingInt;
        %             quadRefine = NumericalQuadrature(vars.refineQuad);
        %             numRefine = quadRefine.num;
        %             if isequal(quadRefine.num,q.num)
        %                 Urefine = this.U;
        %                 Wrefine = this.W;
        %                 Xrefine = X;
        %             else
        %                 [Urefine,~,~,Xrefine] = ...
        %                     dof2point(this.mesh,this.fe_cell,quadRefine.xhat);
        %                 Wrefine = quadRefine.weights(this.mesh);
        %             end
        %             Vh = this;
        %             m = Vh.mesh;
        %             refCell = Vh.cell;
        %             WW = this.W;
        %             nint = size(X,1);
        %             nint_refine = size(Xrefine,1);
        %             ndof = Vh.ndof;
        %             Tab = Vh.dofIndexes;
        %             Nb = refCell.Nb;
        %             [x1,x2,y1,y2] = m.edgesCoords;
        %             UU = this.U;
        %             num = q.num;
        %
        %
        %             X1 = X(:,1);
        %             X2 = X(:,2);
        %             DIFF1 = repmat(X1,1,nint) - repmat(X1',nint,1);
        %             DIFF2 = repmat(X2,1,nint) - repmat(X2',nint,1);
        %             DIFFNorm = sqrt(DIFF1.^2 + DIFF2.^2);
        %             A1 = DIFF1./DIFFNorm.^2;
        %             A2 = DIFF2./DIFFNorm.^2;
        %             A1(or(isnan(A1),isinf(A1))) = 0;
        %             A2(or(isnan(A2),isinf(A2))) = 0;
        %             N = m.normalVector;
        %             N = N./repmat(sqrt(N(:,1).^2 + N(:,2).^2),1,2);
        %             N1 = spdiags(reshape(repmat(N(:,1)',num,1),nint,1),0,nint,nint);
        %             % Replicate the normal vector (first component) for each integration point.
        %             % And make a sparse diag out of it.
        %             N2 = spdiags(reshape(repmat(N(:,2)',num,1),nint,1),0,nint,nint);
        %             Mat = A1*N1+ A2*N2;
        %             Aop = AbstractMatrix(Mat);
        %             %% Correction of singular integrals
        %
        %             singInt = refCell.singularIntegralDictionnary('doubleLayer');
        %
        %             % Correction for the non-derivative part (the k part)
        %
        %             if isnan(threshold)
        %                 % Proximity threshold, under which inegral are considered
        %                 % singular.
        %                 threshold = max(m.length)/2;
        %             end
        %             I = Vh.segmentsCloseToDofs(threshold);
        %
        %             lines = [];
        %             cols = lines;
        %             vals = lines;
        %             lines_refine = [];
        %             cols_refine = lines_refine;
        %             vals_refine = lines_refine;
        %             for dof = 1:ndof
        %                 neighbourSegs = I{dof};
        %                 A1 = repmat((1:num)',1,length(neighbourSegs));
        %                 A2 = repmat((neighbourSegs-1)*num,num,1);
        %                 A1refine = repmat((1:numRefine)',1,length(neighbourSegs));
        %                 A2refine = repmat((neighbourSegs-1)*numRefine,numRefine,1);
        %                 lines_dof = A1(:) + A2(:);
        %                 cols_dof = 0*lines_dof + dof;
        %                 lines_dof_refine = A1refine(:) + A2refine(:);
        %                 cols_dof_refine = 0*lines_dof_refine + dof;
        %                 X_krefine = Xrefine(lines_dof_refine,:);
        %                 vals_dof_refine = zeros(length(lines_dof_refine),Nb);
        %
        %                 for b = 1:Nb;
        %                     seg = find(Tab(:,b)==dof);
        %                     assert(length(seg)<=1);
        %                     if isempty(seg)
        %                         vals_dof_refine(:,b) = 0;
        %                     else
        %                         A = repmat([x1(seg) y1(seg)],length(lines_dof_refine),1);
        %                         B = repmat([x2(seg) y2(seg)],length(lines_dof_refine),1);
        %                         Nloc = repmat(N(seg,:),length(lines_dof_refine),1);
        %                         vals_dof_refine(:,b) = singInt{b}(X_krefine,A, B,Nloc);
        %                     end
        %                 end
        %
        %                 % Old values that need to be substracted :
        %                 l = dof;
        %                 phi_l = UU(:,l);
        %                 p_ind = phi_l~=0;
        %                 W_p = WW(p_ind);
        %                 U_p = phi_l(p_ind);
        %
        %                 Gkp = Mat(lines_dof,p_ind);
        %                 valCorrec = Gkp*(W_p.*U_p);
        %
        %                 % Assembling matrix
        %                 lines = [lines;lines_dof]; % Is it possible to know the
        %                 % size of lines, cols and vals in advance ?
        %                 cols = [cols;cols_dof];
        %                 vals = [vals;-valCorrec];
        %
        %                 lines_refine = [lines_refine;lines_dof_refine];
        %                 cols_refine = [cols_refine;cols_dof_refine];
        %                 vals_refine = [vals_refine;sum(vals_dof_refine,2)];
        %
        %             end
        %             correc = sparse(lines,cols,vals,nint,ndof);
        %             fineQuad = sparse(lines_refine,cols_refine,vals_refine,nint_refine,ndof);
        %
        %             %% Creating Bilinear form
        %             mv_p = @(x)(...
        %                 real(...
        %                 UU'*(WW.*(Aop.mv_prod(WW.*(UU*x)) + correc*x))...
        %                 +	Urefine'*(Wrefine.*(fineQuad*x))));
        %
        %             T = 1/pi*BilinearForm(Vh,mv_p);
        %
        %         end
        % lues that need to be substracted :
        %                 l = dof;
        %                 phi_l = UU(:,l);
        %                 p_ind = phi_l~=0;
        %                 W_p = WW(p_ind);
        %                 U_p = phi_l(p_ind);
        %
        %                 Gkp = Mat(lines_dof,p_ind);
        %                 valCorrec = Gkp*(W_p.*U_p);
        %
        %                 % Assembling matrix
        %                 lines = [lines;lines_dof]; % Is it possible to know the
        %                 % size of lines, cols and vals in advance ?
        %                 cols = [cols;cols_dof];
        %                 vals = [vals;-valCorrec];
        %
        %                 lines_refine = [lines_refine;lines_dof_refine];
        %                 cols_refine = [cols_refine;cols_dof_refine];
        %                 vals_refine = [vals_refine;sum(vals_dof_refine,2)];
        %
        %             end
        %             correc = sparse(lines,cols,vals,nint,ndof);
        %             fineQuad = sparse(lines_refine,cols_refine,vals_refine,nint_refine,ndof);
        %
        %             %% Creating Bilinear form
        %             mv_p = @(x)(...
        %                 real(...
        %                 UU'*(WW.*(Aop.mv_prod(WW.*(UU*x)) + correc*x))...
        %                 +	Urefine'*(Wrefine.*(fineQuad*x))));
        %
        %             T = 1/pi*BilinearForm(Vh,mv_p);
        %
        %         end
        % Deprecated, to be removed
        %         function[R] = hyperSingularOperator(this,varargin)
        %             p = inputParser;
        %             p.addOptional('zeroMeanModif',0)
        %             p.parse(varargin{:});
        %             vars = p.Results;
        %             z = vars.zeroMeanModif;
        %
        %             G = LogKernel(1);
        %             singIntDer = this.fe_cell.singularIntegralDictionnaryDer('ln');
        %             R = -1/(2*pi)*this.BIO('dk',G,'singularIntegralDer',singIntDer,'compressSBD',compressSBD);
        %             if z~= 0
        %                 a = this.meansBasis();
        %                 R = R + z*AbstractMatrix(@(x)(a*(a'*x)),this.ndof,this.ndof);
        %             end
        %         end
        %
        %         function[Mat] = exactCloseIntegrals(this,X,singKernel_id,I,varargin)
        %             % This function computes the matrix
        %             % Aij = \int_[this]K(x_i,y)phi_j(y) if I{j} contains i, 0
        %             % otherwise.
        %             % phi_j is the j-th basis function of the FE space or its
        %             % derivative if the argument 'derivative' is set to true.
        %             % K is the singular kernel represented by singKernel_id.
        %             p = inputParser;
        %             p.addOptional('derivative',false);
        %             p.parse(varargin{:});
        %             vars = p.Results;
        %
        %             nX = size(X,1);
        %             m = this.mesh;
        %             ndof = this.ndof;
        %             T = this.dofIndexes;
        %             feCell = this.fe_cell;
        %             Nb = feCell.Nb;
        %
        %             [Ax,Bx,Ay,By] = m.edgesCoords; % edges [A, B] with A = [Ax, Ay], B = [Bx, By]
        %             if vars.derivative
        %                 singInt = this.fe_cell.singularIntegralDictionnaryDer(singKernel_id);
        %             else
        %                 singInt = this.fe_cell.singularIntegralDictionnary(singKernel_id);
        %             end
        %
        %             lines = []; cols = lines; vals = lines;
        %             dofsToLoop = find(~cellfun('isempty',I))';
        %             % Loop over Y dofs.
        %             for dof_y = dofsToLoop
        %                 lines_X = I{dof_y}';
        %                 cols_dof = 0*lines_X + dof_y;
        %                 X_k = X(lines_X,:);
        %                 val_dofX = zeros(length(lines_X),Nb);
        %                 for b = 1:Nb
        %                     seg = find(T(:,b)==dof_y);
        %                     assert(length(seg)<=1);
        %                     if isempty(seg)
        %                         val_dofX(:,b) = 0;
        %                     else
        %                         AY = repmat([Ax(seg) Ay(seg)],length(lines_X),1);
        %                         BY = repmat([Bx(seg) By(seg)],length(lines_X),1);
        %                         val_dofX(:,b) = singInt{b}(X_k,AY,BY);
        %                     end
        %                 end
        %
        %                 lines = [lines;lines_X];%#ok
        %                 cols = [cols;cols_dof];%#ok
        %                 vals = [vals;sum(val_dofX,2)];%#ok
        %
        %             end
        %             Mat = sparse(lines,cols,vals,nX,ndof);
        %         end
        %         function[Mat] = approxCloseIntegrals(this,X,singKernel_id,I,varargin)
        %             % This function computes the matrix
        %             % Aij = \int_[this]K(x_i,y)phi_j(y) if I{j} contains i, 0
        %             % otherwise, with an approximate integral using the Gauss
        %             % points of this FE space.
        %             % phi_j is the j-th basis function of the FE space or its
        %             % derivative if the argument 'derivative' is set to true.
        %             % K is the singular kernel represented by singKernel_id.
        %             p = inputParser;
        %             p.addOptional('derivative',false);
        %             p.addOptional('constantTerm',false); % Set to true if you need
        %             % to correct only the constant term.
        %             p.parse(varargin{:});
        %             vars = p.Results;
        %             k = singularityDictionary(singKernel_id); % a kernel object.
        %             YGauss = this.gaussPoints;
        %             U = this.phi;
        %
        %             WW = AbstractMatrix.spdiag(this.W);
        %             nX = size(X,1);
        %             ndof = this.ndof;
        %             if vars.derivative
        %                 func = k.der;
        %             else
        %                 func = k.func;
        %             end
        %             lines = []; cols = lines; vals = lines;
        %             dofsToLoop = find(~cellfun('isempty',I))';
        %             % Loop over Y dofs.
        %             for dof_y = dofsToLoop
        %                 lines_X = I{dof_y}';
        %                 cols_dof = 0*lines_X + dof_y;
        %                 X_k = X(lines_X,:);
        %                 l = dof_y;
        %                 phiY_l = U(:,l);
        %                 p_ind = phiY_l~=0;
        %                 WY_p = WW.concretePart(p_ind,p_ind);
        %                 if vars.constantTerm
        %                     UY_p = this.fe_cell.constantTerm(X_k,dof_y); % This is a matrix.
        %                 else
        %
        %                     Y_p_Gauss = YGauss(p_ind,:);
        %                     UY_p = phiY_l(p_ind);
        %                 end
        %                 n_p = nnz(phiY_l);
        %                 n_k = length(lines_X);
        %                 rXkYpGauss = sqrt((repmat(X_k(:,1),1,n_p) - repmat(Y_p_Gauss(:,1)',n_k,1)).^2 ...
        %                     + (repmat(X_k(:,2),1,n_p) - repmat(Y_p_Gauss(:,2)',n_k,1)).^2);
        %                 rXkYpGauss(rXkYpGauss < 1e-8) = 1e-8;
        %                 Gkp = func(rXkYpGauss);
        %                 val_dofX = Gkp*(WY_p*UY_p);
        %
        %                 % Assembling matrix
        %                 lines = [lines;lines_X];%#ok
        %                 cols = [cols;cols_dof];%#ok
        %                 vals = [vals;val_dofX];%#ok
        %
        %             end
        %             Mat = sparse(lines,cols,vals,nX,ndof);
        %         end
        
        
    end
end
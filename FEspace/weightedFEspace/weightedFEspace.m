classdef weightedFEspace < FEspace
    %A space of functions defined on a mesh, endowed with a scalar product
    %of the form {f,g} = \int f(x)g(x)w(x) where w is a weight function.
    %A low order Gaussian quadrature for w is computed on each segment of
    %the mesh.
    
    properties
        I % indexes of the segments where to compute a special quadrature for w.
        weight % the weigth function defined in terms of the curvilinear abscissa
        dweight % derivative of the weight
        weight_id = [];
    end
    
    methods
        function[w,dw] = FEweight(Vh,id)
            switch id
                case 'sqrt(1-t^2)'
                    L = sum(Vh.mesh.length);
                    w = @(s)(sqrt(s.*(L - s)));
                    dw = @(s)((L/2-s)./w(s));
                case '1/sqrt(1-t^2)'
                    L = sum(Vh.mesh.length);
                    w = @(s)(1./sqrt(s.*(L-s)));
                    dw = @(s)((s - L/2).*w(s).^3);
                otherwise
                    error('unknown weight function');
                    
            end
        end
        function[this] = weightedFEspace(m,c,w,varargin)
            % Parse inputs
            p = inputParser;
            p.addOptional('quadNum',3);
            p.addOptional('specialQuadSegs',[1:5 (m.nseg-4:-1:m.nseg)]);
            p.parse(varargin{:});
            % Instantiate usinge FEspace class.
            q = p.Results.quadNum; % number of quadrature points on each seg.
            this@FEspace(m,c,q);
            
            % Define the weight
            this.weight_id =w; % a string
            [this.weight,this.dweight] = this.FEweight(w);
            w = this.weight; % a handle function
            
            % Update the quadrature weights to account for the weight
            sVec = this.sVec;
            sVert = m.sVertices;
            this.W = this.W.*w(sVec);
            
            % Now compute a new quadrature for semgents near the
            % singularities of the weight. We use the variable change
            % s = L/2*(1 + cos(theta))
            T = this.dofIndexes;
            l = m.length;
            L = sum(l);
            I = p.Results.specialQuadSegs;
            % elementary quadrature after variable change
            [uhat,u_what] = Gauss_Legendre1D(q,0,1);
            uhat = flipud(uhat); u_what = flipud(u_what);
%             for i = 1:length(I) % loop over the segments where the quad
%                 % needs to be changed.
%                 seg_num = I(i);
%                 AB = m.segment(seg_num);
%                 A = AB(1); B = AB(2); % indexes of the vertex defining A and B
%                 sA = sVert(A);
%                 sB = sVert(B);
%                 % get coords of A and B.
%                 A = this.mesh.vertex(A);
%                 B = this.mesh.vertex(B);
%                 % transformation  s = L/2*(1 + cos(theta))
%                 % We let x = cos(theta), s = L/2(1 + x)
%                 theta1 = real(acos(2*sA/L-1));
%                 theta2 = real(acos(2*sB/L-1));
%                 xhat = sort(cos(theta1 + (theta2 - theta1)*uhat));
%                 if strcmp(this.weight_id,'1/sqrt(1-t^2)')
%                     what = abs((theta2 - theta1))*u_what;
%                 else
%                     what = (1-xhat.^2).*abs((theta2 - theta1)).*u_what;
%                 end
%                 sref = ((xhat +1)*L/2 - sA)/(sB-sA); % map to the unit segment.
%                 idx = (seg_num-1)*q + (1:q)';
%                 u = (B-A);
%                 this.gaussPoints(idx,1) = A(1) + u(1)*sref;
%                 this.gaussPoints(idx,2) = A(2) + u(2)*sref;
%                 for b = 1:this.fe_cell.Nb
%                     [phi,dxphi] = this.fe_cell.basis(sref);
%                    
%                     dof = T(seg_num,b);
%                     this.phi(idx,dof) = phi(:,b);
%                     this.dphi(idx,dof) = dxphi(:,b)/l(seg_num);
%                 end
%                 this.W(idx) = what;
%                 this.sVec(idx) = sA + sref*(sB - sA);
%             end
            % Vectorial version of the loop:
            T = this.dofIndexes; 
            sA = sVert(T(:,1));
            sB = sVert(T(:,2));
            A = this.mesh.vertices(T(:,1),:);
            B = this.mesh.vertices(T(:,2),:);
            
            theta1 = real(acos(2*sA/L-1));
            theta2 = real(acos(2*sB/L-1));
            xhat = reshape(cos(theta1 + (theta2 - theta1)*uhat')',this.nint,1);
            if strcmp(this.weight_id,'1/sqrt(1-t^2)')
                what = repelem(abs(theta2 - theta1),this.quad.num).*repmat(u_what,this.nseg,1);
            else
                what = (1-xhat.^2).*repelem(abs(theta2 - theta1),this.quad.num).*repmat(u_what,this.nseg,1);
            end
            sA = repmat(sA,1,this.quad.num)';
            sA = sA(:);
            sB = repmat(sB,1,this.quad.num)';
            sB = sB(:);
            sref = ((xhat +1)*L/2 - sA)./(sB-sA); % map to the unit segment.
            this.sVec = sA + (sB-sA).*sref;
            u = B - A;
            u = repelem(u,this.quad.num,1);
            A = repelem(A,this.quad.num,1);
            gP1 = A(:,1) + u(:,1).*sref;
            gP2 = A(:,2) + u(:,2).*sref;
            this.gaussPoints = [gP1 gP2];
            this.phi = 0*this.phi;
            this.dphi = 0*this.dphi;
            [phi,dxphi] = this.fe_cell.basis(sref);
            for b = 1:this.fe_cell.Nb
                lines = repmat((1:this.quad.num),1,this.nseg);
                lines = lines + repelem((0:this.nseg-1),this.quad.num)*this.quad.num;
                cols = repelem(T(:,b),this.quad.num);
                this.phi = this.phi + sparse(lines,cols,phi(:,b),this.nint,this.ndof);
                this.dphi = this.dphi + sparse(lines,cols,dxphi(:,b)./repelem(l,this.quad.num),this.nint,this.ndof);
            end
            this.W = what;
            Wdiag = AbstractMatrix.spdiag(this.W);
            M_mat = AbstractMatrix(this.phi'*Wdiag*this.phi);
            dM_mat = AbstractMatrix(this.dphi'*Wdiag*this.dphi);
            this.Mass = BilinearForm(this,this,M_mat);
            this.dMass = BilinearForm(this,this,dM_mat);
            this.selectedTest = this.phi;
            this.I = I;
            this.weight = w;
        end
        function[this] = remesh(this,N)
            m = this.mesh.remesh(N);
            if ~isempty(this.weight_id)
                this = weightedFEspace(m,this.fe_cell.id,this.weight_id,this.quad.num);
            else
                this = weightedFEspace(m,this.fe_cell.id,this.weight,this.quad.num);
            end
            
        end
        function[] = showWeight(this)
            s = this.sVec;
            plot(s,this.weight(s));
            xlabel('curvilinear absissa')
            title(sprintf('Weight function %s',this.weight_id));
        end
        function[Mat] = omega_dx_omega(this)
            % Returns the matrix that maps Aij such that
            % Aij = (\omega dx \omega phi_i)(x_j) where phi_i is the i-th
            % basis function and x_j is the j-th gauss point.
            L = sum(this.mesh.length);
            assert(strcmp(this.weight_id,'1/sqrt(1-t^2)'));
            omega2 = 1./(this.weight(this.sVec)).^2;
            omega2(or(isinf(omega2),isnan(omega2))) = 0;
            s = this.sVec - L/2; % because dw(s) = (L/2 - s)/w(s)
            omega2_diag = spdiags(omega2,0,this.nint,this.nint);
            s_diag =  spdiags(s,0,this.nint,this.nint);
            Mat = omega2_diag*this.dphi - s_diag*this.phi;
        end
        function[Mat] = omega2(this)
            if strcmp(this.weight_id,'sqrt(1-t^2)')
                om2 = this.weight(this.sVec).^2;
                om2(or(isinf(om2),isnan(om2))) = 0;
                om2_diag = spdiags(om2,0,this.nint,this.nint);
                Wd = spdiags(this.W,0,this.nint,this.nint);
                Mat = this.phi'*Wd*om2_diag*this.phi;
            else
                om2 = 1./this.weight(this.sVec).^2;
                om2(or(isinf(om2),isnan(om2))) = 0;
                om2_diag = spdiags(om2,0,this.nint,this.nint);
                Mat = om2_diag*this.phi;
            end
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
            switch p.Results.correcMethod
                case 'constantTerm'
                    Mat = this.correctionConstantTerm(X,singKernel_id);
                case 'precise'
                    Mat  = this.correctionPrecise(X,singKernel_id,varargin);
            end
        end
        function[Mat] = ln_omega_reg(this,X)
            % Vh a weighted FE space.
            % This function computes the sparse matrix
            % Aik = Cik(\int_[this]K(X_k,y)  - \intApprox_[this]K(X_k,y))
            % where the non-zeros indexes are chosen when X_k is close to
            % the segments attached to the dof number i. Cik is the
            % piecewise constant function on each segment of the mesh such
            % that the i-th basis function can be written on each segment
            % as phi_i(Y) = Cik(Y) + someConstant times (X_k - Y)
            
            M = size(X,1); % Number of points in X
            me = this.mesh;
            l = me.length; % vector containing length of each segment in the mesh.
            L = sum(l);
            ndof = this.ndof;
            T = this.dofIndexes;
            feCell = this.fe_cell;
            Nb = feCell.Nb;
            singK = logSingK;
            func = singK.k.func;
            coords = this.mesh.edgesCoords;
            Ax = coords(:,1); Ay = coords(:,2); Bx = coords(:,3); By = coords(:,4);
%             [Ax,Bx,Ay,By] = this.mesh.edgesCoords; % edges [A, B] with A = [Ax, Ay], B = [Bx, By]
            lines = cell(Nb,1); cols = cell(Nb,1); vals = cell(Nb,1);
            %             threshold = 4*max(l);
            %             I = isClose(X,this.dofCoords,threshold);
            %             dofsToLoop = find(~cellfun('isempty',I))';
            for i = 1:this.ndof %dofsToLoop
                for b = 1:Nb
                    segNum = find(T(:,b)==i); % On which segment is
                    % the b-th basis function attached to dof i.
                    %                     X_k = X(I{i}',:);
                    if isempty(segNum)
                        % nothing to do
                    else
                        A = [Ax(segNum) Ay(segNum)]; B = [Bx(segNum) By(segNum)]; % the segment is [A,B].
                        
                        
                        Mid = (A + B)/2; % segment middle.
                        lAB2 = l(segNum)^2*4; % squared length of the segment.
                        
                        
                        
                        % We select in X the points that are within three balls centered
                        % respectively on A,B and Mid, and of radius lAB.
                        AX2 = (X(:,1) - Ax(segNum)).^2 + (X(:,2) - Ay(segNum)).^2;
                        BX2 = (X(:,1) - Bx(segNum)).^2 + (X(:,2) - By(segNum)).^2;
                        MX2 = (X(:,1) - Mid(1)).^2 + (X(:,2) - Mid(2)).^2;
                        ks = find(or(AX2 < lAB2, or(BX2 < lAB2,MX2<lAB2)));
                        X_k = X(ks,:); % seclected points in X.
                        %for each k, parameter cb(k) such that phi_b(Y) = cb(k) + C*(X_k(k,:) - Y)
                        % where C is some constant.
                        approxInt = this.I0approx(func,X_k,A,B,segNum);
                        cb = this.fe_cell.constantTerm(b,X_k,A,B);
                        db = this.fe_cell.linearTerm(b,X_k,A,B);
                        [~,~,~,~,~,~,alpha,~,~] = parameters_singInt(X_k,repmat(A,size(X_k,1),1),repmat(B,size(X_k,1),1));
                        [sA,~] = this.mesh.s_seg(segNum);
                        s0 = sA - alpha;
                        omega2 = 1./this.weight(s0).^2;
                        omega2(or(isinf(omega2),isnan(omega2))) = 0;
                        H = (L/2-s0).*cb + omega2.*db;
                        % Compute the approximation that was used for
                        % \int_{[A,B]}G(X_k(k,:),Y) dY for each k.
                        
                        
                        % Compute the exact integral to replace :
                        exactInt = this.I0exact(singK, X_k,A,B,segNum); % method that computes
                        % exactly \int_{[A,B]} G(X_k(k,:),Y) dY for each k.
                        
                        % Store the correction.
                        correc_kib = H.*(exactInt - approxInt);
                        lines{b} = [lines{b};ks];
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
        function[Mat] = ln_omega2_reg(this,X)
            % Vh a weighted FE space.
            % This function computes the sparse matrix
            % Aik = Cik(\int_[this]K(X_k,y)  - \intApprox_[this]K(X_k,y))
            % where the non-zeros indexes are chosen when X_k is close to
            % the segments attached to the dof number i. Cik is the
            % piecewise constant function on each segment of the mesh such
            % that the i-th basis function can be written on each segment
            % as phi_i(Y) = Cik(Y) + someConstant times (X_k - Y)
            
            M = size(X,1); % Number of points in X
            N = this.mesh.normalVector;
            me = this.mesh;
            l = me.length; % vector containing length of each segment in the mesh.
%             L = sum(l);
            ndof = this.ndof;
            T = this.dofIndexes;
            feCell = this.fe_cell;
            Nb = feCell.Nb;
            singK = logSingK;
            func = singK.k.func;
            
            EDG = this.mesh.edgesCoords; % edges [A, B] with A = [Ax, Ay], B = [Bx, By]
            Ax = EDG(:,1); Ay = EDG(:,2); Bx = EDG(:,3); By = EDG(:,4);
            lines = cell(Nb,1); cols = cell(Nb,1); vals1 = cell(Nb,1);vals2 = cell(Nb,1);
            %             threshold = 4*max(l);
            %             I = isClose(X,this.dofCoords,threshold);
            %             dofsToLoop = find(~cellfun('isempty',I))';
            for i = 1:this.ndof %dofsToLoop
                for b = 1:Nb
                    segNum = find(T(:,b)==i); % On which segment is
                    % the b-th basis function attached to dof i.
                    %                     X_k = X(I{i}',:);
                    if isempty(segNum)
                        % nothing to do
                    else
                        
                        A = [Ax(segNum) Ay(segNum)]; B = [Bx(segNum) By(segNum)]; % the segment is [A,B].
                        NAB = N(segNum,:);
                        NAB = NAB./(norm(NAB));
                        
                        Mid = (A + B)/2; % segment middle.
                        lAB2 = l(segNum)^2*2; % squared length of the segment.
                        
                        
                        
                        % We select in X the points that are within three balls centered
                        % respectively on A,B and Mid, and of radius lAB.
                        AX2 = (X(:,1) - Ax(segNum)).^2 + (X(:,2) - Ay(segNum)).^2;
                        BX2 = (X(:,1) - Bx(segNum)).^2 + (X(:,2) - By(segNum)).^2;
                        MX2 = (X(:,1) - Mid(1)).^2 + (X(:,2) - Mid(2)).^2;
                        ks = find(or(AX2 < lAB2, or(BX2 < lAB2,MX2<lAB2)));
                        X_k = X(ks,:); % seclected points in X.
                        %for each k, parameter cb(k) such that phi_b(Y) = cb(k) + C*(X_k(k,:) - Y)
                        % where C is some constant.
                        [~,~,~,~,~,~,alpha,~,~] = parameters_singInt(X_k,repmat(A,size(X_k,1),1),repmat(B,size(X_k,1),1));
                        [sA,~] = this.mesh.s_seg(segNum);
                        s0 = sA - alpha;
                        approxInt = this.I0approx(func,X_k,A,B,segNum);
                        cb = this.fe_cell.constantTerm(b,X_k,A,B);
                        
                        
                        
                        omega2 = 1./this.weight(s0).^2;
                        omega2(or(isinf(omega2),isnan(omega2))) = 0;
                        H = omega2.*cb;
                        % Compute the approximation that was used for
                        % \int_{[A,B]}G(X_k(k,:),Y) dY for each k.
                        
                        
                        % Compute the exact integral to replace :
                        exactInt = this.I0exact(singK, X_k,A,B,segNum); % method that computes
                        % exactly \int_{[A,B]} G(X_k(k,:),Y) dY for each k.
                        
                        % Store the correction.
                        correc_kib = H.*(exactInt - approxInt);
                        lines{b} = [lines{b};ks];
                        cols{b} = [cols{b};0*ks + i];
                        vals1{b} = [vals1{b};correc_kib*NAB(1)];
                        vals2{b} = [vals2{b};correc_kib*NAB(2)];
                    end
                   
                    
                end
            end
            % We have a corrective patch for each b on the whole mesh. Sum
            % them to get the final correction.
            Mat1 = sparse(lines{1},cols{1},vals1{1},M,ndof);
            Mat2 = sparse(lines{1},cols{1},vals2{1},M,ndof);
            for b = 2:Nb
                Mat1 = Mat1 + sparse(lines{b},cols{b},vals1{b},M,ndof);
                Mat2 = Mat2 + sparse(lines{b},cols{b},vals2{b},M,ndof);
            end
            Mat = cell(2,1);
            Mat{1} = Mat1;
            Mat{2} = Mat2;
        end
        function[res] = I0approx(this,func,X_k,A,B,segNum)
            [~,~,~,~,~,~,alpha,~,~] = parameters_singInt(X_k,repmat(A,size(X_k,1),1),repmat(B,size(X_k,1),1));
            [sA,~] = this.mesh.s_seg(segNum);
            s0 = sA - alpha;
            q = this.quad.num;
            L = sum(this.mesh.length);
            % Gauss points and weights on segment
            idx = (segNum-1)*q + (1:q)';
            sq = this.sVec(idx);
            Wq = this.W(idx); % weights
            phi = @(s)(-real(acos(2*s/L - 1)));
            Akq = (phi(sq)' - phi(s0))./(repmat(this.weight(s0),1,length(sq)));
            Akq(abs(Akq)<1e-15) = 1e-15;
            Gkq = func(abs(Akq));
            
            res = (Gkq*Wq);
        end
        function[res] = I0exact(this,singK,X_k,A,B,segNum)
            F = singK.primitiveOfKernel;
            % Gauss points and weights on segment
            L = sum(this.mesh.length);
            [~,~,~,~,~,~,alpha,~,~] = parameters_singInt(X_k,repmat(A,size(X_k,1),1),repmat(B,size(X_k,1),1));
            [sA,sB] = this.mesh.s_seg(segNum);
            s0 = sA - alpha;
            phi = @(s)(-real(acos(2*s/L - 1)));
            res = this.weight(s0).*(F((phi(sB) - phi(s0))./(this.weight(s0)))...
                - F((phi(sA) - phi(s0))./(this.weight(s0))));
        end
        function[Mat] = correctionPrecise(this,X,singKernel_id,I,varargin)
            % Vh a weighted FE space.
            % This function computes the matrix
            % Aij = \int_[this]K(x_i,y)phi_j(y) if I{j} contains i, 0
            % otherwise.
            % phi_j is the j-th basis function of the FE space or its
            % derivative if the argument 'derivative' is set to true.
            % singKernel is the singular kernel. It must be able to return the handle
            % function corresponding to the singular kernel, and moreover must also
            % provide a function I(X,A,B) that computes \int_{[A,B]}\ln(|X - Y|)dY
            % for any list of M R   A2tildex = w0.*(Gkq*(Wq./this.weight(sq)));^2 points (X1,...XM) and for two R^2 points A and B
            
            M = size(X,1);
            me = this.mesh;
            L = sum(me.length);
            ndof = this.ndof;
            T = this.dofIndexes;
            feCell = this.fe_cell;
            Nb = feCell.Nb;
            q = this.quad.num;
            gaussY = this.gaussPoints;
            WY = this.W;
            singK = singKernel(singKernel_id);
            func = singK.k.func;
            [Ax,Bx,Ay,By] = me.edgesCoords; % edges [A, B] with A = [Ax, Ay], B = [Bx, By]
            lines = []; cols = lines; vals = lines;
            dofsToLoop = find(~cellfun('isempty',I))';
            % Loop over Y dofs.
            for i = dofsToLoop
                lines_X = I{i}';
                
                cols_dof = 0*lines_X + i;
                X_k = X(lines_X,:);
                Mk = size(X_k,1);
                correc_ki = zeros(length(lines_X),Nb);
                for b = 1:Nb
                    seg = find(T(:,b)==i);
                    assert(length(seg)<=1);
                    if isempty(seg)
                        correc_ki(:,b) = 0;
                    else
                        A = [Ax(seg) Ay(seg)];
                        B = [Bx(seg) By(seg)];
                        cb = this.fe_cell.constantTerm(b,X_k,A,B);
                        db = this.fe_cell.linearTerm(b,X_k,A,B);
                        [~,~,~,~,~,u,alpha,beta,d] = parameters_singInt(X_k,repmat(A,size(X_k,1),1),repmat(B,size(X_k,1),1));
                        [sA,sB] = this.mesh.s_seg(seg);
                        s0 = sA -alpha;
                        s0(s0 < 0) = -s0(s0 < 0);
                        s0(s0>L) = 2*L-s0(s0 > L);
                        idx = (seg-1)*q + (1:q)';
                        Yq = gaussY(idx,:);
                        Wq = WY(idx);
                        w0 = this.weight(sA);
                        dw0 = this.dweight(sA);
                        rXkYq = sqrt((repmat(X_k(:,1),1,q) - repmat(Yq(:,1)',Mk,1)).^2 ...
                            + (repmat(X_k(:,2),1,q) - repmat(Yq(:,2)',Mk,1)).^2);
                        rXkYq(rXkYq < 1e-15) = 1e-15;
                        Gkq = func(rXkYq);
                        %                         Gkq(or(isnan(Gkq),isinf(Gkq))) = 0;
                        sq = this.sVec(idx);
                        sq_s0 = repmat(sq',length(s0),1)-repmat(s0,1,q);
                        if ~ismember(seg,this.I)
                            % in this case, the only singularity comes
                            % from the kernel,
                            A2x_constant = cb.*w0.*(singK.I0seg(X_k,A,B));
                            A2tildex_constant = cb.*w0.*(Gkq*(Wq./this.weight(sq)));
                            correc_constant = A2x_constant - A2tildex_constant;
                            A2x_linear = (db.*w0 + cb*dw0).*(1/2*singK.I1(alpha,beta,d));
                            
                            A2tildex_linear = (db.*w0 + cb*dw0).*((Gkq.*sq_s0)*(Wq./this.weight(sq)));
                            correc_linear = A2x_linear - A2tildex_linear;
                            correc_ki(:,b) = correc_constant + correc_linear;
                        else
                            % In this case, we have to do things manually.
                            A2tildex_constant = cb.*(Gkq*Wq);
                            A2tildex_linear = db.*((Gkq.*sq_s0)*Wq);
                            A2x_constant = 0*A2tildex_constant;
                            A2x_linear = 0*A2tildex_linear;
                            for kk = 1:Mk
                                x = X_k(kk,:);
                                y1 = @(ss)(A(1) + (ss-sA)*u(kk,1));
                                y2 = @(ss)(A(2) + (ss-sA)*u(kk,2));
                                Y = @(ss)([y1(ss(:)) y2(ss(:))]);
                                r = @(ss)(sqrt(cWise_dot(repmat(x,length(ss),1) - Y(ss(:)),...
                                    repmat(x,length(ss),1) - Y(ss(:))))');
                                dr = @(ss)(rprime(r,x,Y,u(kk,:),ss(:)));
                                fun = @(theta)((r((cos(theta) + 1)*L/2)));
                                dfun = @(theta)((-dr((cos(theta) + 1)*L/2).*sin(theta)*L/2));
                                theta_a = real(acos(2*sA/L - 1));
                                theta_b = real(acos(2*sB/L - 1));
                                theta_0 = real(acos(2*s0(kk)/L - 1));
                                phi = @(ss)(cos(ss(:))' - cos(theta_0));
                                phi_der = [-sin(theta_0);-cos(theta_0);sin(theta_0)];
                                if dfun(theta_0) ~= 0
                                    A2x_constant(kk,1) = cb(kk)*singK.Icomp_r(fun,theta_0,dfun(theta_0),theta_b,theta_a);
                                    A2x_linear(kk,1) = db(kk)*singK.Icomp_r_phi(fun,theta_0,dfun(theta_0),theta_b,theta_a,phi,phi_der);
                                else
                                    cou
                                end
                                
                                
                            end
                            correc_ki(:,b) = A2x_constant - A2tildex_constant + A2x_linear - A2tildex_linear;
                        end
                    end
                    
                end
                
                lines = [lines;lines_X];%#ok
                cols = [cols;cols_dof];%#ok
                vals = [vals;sum(correc_ki,2)];%#ok
                
            end
            Mat = sparse(lines,cols,vals,M,ndof);
            function[rp] = rprime(r,x,Y,u,s)
                x = repmat(x,length(s),1);
                u = repmat(u,length(s),1);
                rp = cWise_dot(Y(s(:))-x,u)./r(s(:));
                rp(or(isnan(rp),isinf(rp))) = 1;
            end
        end
    end
end



classdef Op < AbstractMatrix
    % Op objects represent a matrix A which elements are given by
    % A(i,j) = G(|X(i) - Y(j)|) for some set of points X and Y. However, this
    % object does not store the full matrix, instead it uses a Fourier
    % approximation to reduce the amount of data stored and the time needed
    % for computing a Matrix-Vector product.
    properties
        rMax = 0; % when the user inputs X and Y as data set, rMax is an upper
        % bound of the maximal distance between a point in X and a point in
        % Y. The data set are then rescaled by X/rMax and Y/rMax
        x = [];
        y = []; % rescaled data
        kernel@Kernel = Kernel;% Kernel object. It contains the function G and other information needed for the method
        a = 1; % Bounds of the radial quadrature
        tol = 0; % required accuracy of the approximation
        rq@RadialQuadrature = RadialQuadrature;
        q2d@Quad2D = Quad2D; % Quad2D object, computes the far interaction
        timesClose = struct('total',0,'computeInteractions',0,'computingTree',NaN,'rangeSearch',NaN,...
            'NUFFTcloseField',NaN,'assemble',0);
        timeTotalAssembling = 0; % Timing info
        full = false;
        Y_kdtree = [];
        rxy;
        verbose;
    end
    
    methods (Access = public)
        %% Constructor
        function[op] = Op(X,k,Y,varargin)
            op@AbstractMatrix([],[],0,0);
            tTotalAssembling = tic;
            if nargin == 0
                % Empty Op.
                return
            end
            %% Argument Parsing
            p = inputParser;
            p.KeepUnmatched = true;
            p.addOptional('noCloseField',false); % Set to true to disable close field
            % This is useful when you know that clouds XX and YY are far
            % apart.
            p.addOptional('principalPart',false); % Turn this on if you
            % want to compute the Operator
            % A(i,jà = G_loc|X(i) - Y(j)|
            % Where G_loc is such that there exists a C^\infty function R
            % with
            % G_loc + R = G
            p.addOptional('fullMatrix',false);
            p.addOptional('a',Op.default_a(size(X,1),size(Y,1))); p.addOptional('a_factor',1);
            p.addOptional('tol',1e-8); p.addOptional('precomputedRadialQuad',[]);
            p.addOptional('precomputedKDtree',[]);
            p.addOptional('rMax',Op.rMaxCalc(X,Y));
            p.addOptional('verbose',1);
            % If user wishes to compute a radial quad with specific option,
            % or use several times the same radialquad for different Op
            % object, this is the way to pass it as an argument.
            p.parse(varargin{:});
            optVars = p.Results;
            op.Y_kdtree = optVars.precomputedKDtree;
            
            %% Preparing data and rescaling
            op.rMax = optVars.rMax;
            op.x = X/op.rMax;
            op.y = Y/op.rMax;
            op.N1 = size(X,1);
            op.N2 = size(Y,1);
            op.kernel = k.dilatation(op.rMax);
            op.verbose = optVars.verbose;
            %% Create the operator
            if optVars.fullMatrix
                % The operator is fully assembled.
                
                % the timing information are normally set by the call to
                % Radial
                op = op.computeFullMatrix;
                % Computes all the Nx*Ny interactions
                % Also sets the timesClose struct.
            else
                % A radial quadrature is performed on the kernel for the
                % far interactions
                % Local interactions are computed explicitly
                op.tol = optVars.tol;
                op.a = optVars.a_factor*optVars.a;
                
                %% Far operator
                if isempty(optVars.precomputedRadialQuad)
                    op.rq = op.kernel.radialQuadKernel(op.a,op.tol,varargin{:});
                else
                    op.rq = optVars.precomputedRadialQuad;
                end
                op.q2d = Quad2D(op.rq);
                
                %% Local Part
                if ~or(isa(op.kernel,'J0Kernel'),optVars.noCloseField)
                    op = createLocalPart(op);
                    % Computes the local interactions
                    % Also sets the timesClose struct;
                else
                    % No need to compute the close interactions.
                    op.concretePart = sparse([],[],[],op.N1,op.N2);
                    if optVars.noCloseField
                        rMinTest = op.rMinTest;
                        if rMinTest < op.a*op.rMax
                            warning(['There are pairs of points closer than ', ...
                                'parameter ''a''. It was however requested not ',...
                                'to compute the local interactions. This may lead to',...
                                ' inaccurate results']);
                        end
                    end
                end
                if optVars.principalPart
                    % Erase the radial quadratue components
                    op.rq = RadialQuadrature;
                    op.q2d = Quad2D;
                end
                op.abstractPart = @(vec)(op.q2d.conv(op.x,op.y,vec));
            end
            
            op.timeTotalAssembling = toc(tTotalAssembling);
        end
        %% Immediate properties
        function[nxi] = Nxi(this,test_a)
            if nargin == 1
                nxi = this.q2d.Nxi;
            else
                rq_test = this.kernel.radialQuadKernel(test_a,this.tol);
                nxi = sum(fix(exp(1)/2*0.8*rq_test.rho + 5*log(1/(this.tol+10^(-10)))));
            end
        end
        function[out] = NCI(op,test_a)
            if nargin == 1
                out = nnz(op.concretePart);
            else
                if isempty(op.Y_kdtree)
                    I = rangesearch(op.Y,op.X,op.rMax*test_a*1.05);
                else
                    I = rangesearch(op.Y_kdtree,op.X,op.rMax*test_a*1.05);
                end
                out = length(cell2mat(I')');
            end
        end
        function[tFar] = timesFar(this)
            tFar.radialQuadrature = this.rq.times.total;
            tFar.quad2D = this.q2d.time;
            tFar.total = tFar.radialQuadrature + tFar.quad2D;
        end
        function[r] = rMinTest(this)
            % First pass
            x1 = mean(this.x,1);
            [~,ind] = min(sqrt((x1(1) - this.y(:,1)).^2 + (x1(2) - this.y(:,2)).^2));
            y1 = this.y(ind,:);
            % Second pass
            [~,ind] = min(sqrt((this.x(:,1) - y1(1)).^2 + (this.x(:,2) - y1(2)).^2));
            x2 = this.x(ind,:);
            [~,ind] = min(sqrt((x2(1) - this.y(:,1)).^2 + (x2(2) - this.y(:,2)).^2));
            y2 = this.y(ind,:);
            % Third pass
            [~,ind] = min(sqrt((this.x(:,1) - y2(1)).^2 + (this.x(:,2) - y2(2)).^2));
            x3 = this.x(ind,:);
            [~,ind] = min(sqrt((x3(1) - this.y(:,1)).^2 + (x3(2) - this.y(:,2)).^2));
            y3 = this.y(ind,:);
            r = norm(x3 - y3,2)*this.rMax;
        end
        function[X] = X(op)
            X = op.x * op.rMax;
        end
        function[Y] = Y(op)
            Y = op.y * op.rMax;
        end
        function[out] = a_factor(op)
            % Size of parameter a compared to default value.
            out = op.a/Op.default_a(op.N1,op.N2);
        end
        
        %% Setters
        function[op] = createLocalPart(op)
            times.rangeSearch = 0;
            times.computeInteractions = 0;
            times.NUFFTcloseField = 0;
            times.assemble = 0;
            tTotalTime = tic;
            % Find close interactions
            tRangeSearch = tic;
            times.computingTree = 0;
            if isempty(op.Y_kdtree)
                tComputingTree = tic;
                op.Y_kdtree = KDTreeSearcher(op.Y);
                times.computingTree = toc(tComputingTree);
            end
            [I,rxyTemp] = rangesearch(op.Y_kdtree,op.X,op.rMax*op.a*1.05);
            jdx = cell2mat(I')';
            op.rxy = cell2mat(rxyTemp')'/op.rMax;
            idx = zeros(size(jdx));
            sp_ind = 1;
            for x_ind=1:length(I)
                idx(sp_ind:(sp_ind+length(I{x_ind})-1)) = x_ind;
                sp_ind = sp_ind + length(I{x_ind});
            end
            % Save memory
            times.rangeSearch = toc(tRangeSearch);
            NCI = length(op.rxy);
%             Nxi = op.q2d.Nxi;
            % close-field
            if NCI ~= 0
                rxyApply = op.rxy;
                rxyApply(rxyApply*op.rMax < 1e-15) = 1e-15/op.rMax; 
                exact_Interactions = op.kernel.func(rxyApply); % Exact local interactions
                
%                 exact_Interactions(or(isinf(exact_Interactions),isnan(exact_Interactions))) = 0;
                tAssemble = tic;
                C1 = sum(abs(op.rq.alpha0.*Cp(op.rq.rho)).*op.rq.rho.^2); % Bound for the second derivative.
                Ninterp = fix(sqrt(C1)*op.a/sqrt(8*op.tol))+10; % Guarantees interpolation error < tol.
                xinterp = linspace(0,op.a*1.1,Ninterp);
                yinterp = op.rq.eval(xinterp);
                radial_quadratureNear0 = interp1(xinterp,yinterp,op.rxy); % we remove the radial
                % quadrature contribution (using interpolation to avoid evaluating at
                % all points).
                
                C_val = exact_Interactions - radial_quadratureNear0;
                
                op.concretePart = sparse(idx,jdx,C_val,op.N1,op.N2);
%                 % local correction of error due to NUFFT
%                 if op.verbose >=1
%                     fprintf('NUFFT for local correction : *\n')
%                 end
%                 tNUFFTcloseField = tic;
%                 B1_inds = op.kernel.func(rxyApply);   
% %                 B1_inds(or(isinf(B1_inds),isnan(B1_inds))) = 0;
%                 B2_inds = nufft2d3(Nxi, op.q2d.xi_nu(:,1), op.q2d.xi_nu(:,2), ...
%                     op.q2d.w_nu, +1, op.tol/100, length(idx),...
%                     op.x(idx,1) - op.y(jdx,1), op.x(idx,2) - op.y(jdx,2));
%                 times.NUFFTcloseField = toc(tNUFFTcloseField);
%                 % Sparse matrix
%                 tAssemble = tic;
%                 B_inds = B1_inds - B2_inds - op.q2d.offset;
%                 op.concretePart = sparse(idx,jdx,B_inds,op.N1,op.N2);
                times.assemble = toc(tAssemble);
            else
                % No close interactions
                op.concretePart = sparse(op.N1,op.N2); % all zeros
            end
            
            times.total = toc(tTotalTime);
            op.timesClose = times;
            
        end
        function[op] = update_X(op,X,varargin)
            tUpdating = tic;
            p = inputParser;
            p.addOptional('a',op.a);
            p.addOptional('tol',op.tol);
            p.addOptional('principalPart',false);
            p.addOptional('a_factor',1);
            p.parse(varargin{:});
            % Replace the operator by A(i,j) = G(|X(i) - Y(j)|) with new X
            % data.
            % 1°) New rescaling
            rMaxNew = Op.rMaxCalc(X,op.Y);
            if any([rMaxNew > op.rMax,op.tol > p.Results.tol, p.Results.a*p.Results.a_factor ~= op.a,1])
                % Then we have to update the quadrature since it is not
                % valid up to rMaxTest or to requested tolerance
                disp('Recomputing a quadrature');
                rMaxPrevious = op.rMax;
                op.rMax = rMaxNew;
                op.x = X/op.rMax;
                op.y = op.y*rMaxPrevious/op.rMax;
                op.a = max(min(p.Results.a*p.Results.a_factor,0.5),0.9*op.rMinTest/op.rMax);
%                 op.a = p.Results.a*p.Results.a_factor;
                op.N1 = size(X,1);
                op.kernel = op.kernel.dilatation(rMaxNew/rMaxPrevious);
                if op.full
                    op = op.computeFullMatrix;
                else
                    op.a = max(min(p.Results.a*p.Results.a_factor,0.5),0.9*op.rMinTest/op.rMax);
                    op.rq = op.kernel.radialQuadKernel(op.a,op.tol);
                    op.q2d = Quad2D(op.rq);
                    if p.Results.principalPart
                        % Erase the radial quadratue components
                        op.rq = RadialQuadrature;
                        op.q2d = Quad2D;
                    end
                end
            else
                % We can reuse the previous quadrature
                op.x = X/op.rMax;
                op.N1 = size(X,1);
            end
            op = op.createLocalPart;
            op.abstractPart = @(vec)(op.q2d.conv(op.x,op.y,vec));
            op.timeTotalAssembling = toc(tUpdating);
        end
        function[op] = change_a(op,new_a)
            op.a = new_a;
            op.rq = op.kernel.radialQuadKernel(op.a,op.tol);
            op.q2d = Quad2D(op.rq);
            %% Local Part
            % If we decrease a, we just have to remove the pairs that are a
            % bit too far.
            if ~isa(op.kernel,'J0Kernel')
                op = op.createLocalPart;
                % Computes the local interactions
                % Also sets the timesClose struct
            end
            op.abstractPart = @(vec)(op.q2d.conv(op.x,op.y,vec));
        end
        function[op] = balanceMemory(op)
            Nxi = op.Nxi;
            NCI = op.NCI;
            NCI0 = NCI;
            memClose = NCI*16;
            memFar = Nxi*24;
            alpha = 1; % Guessed dependance of NCI in variable a.
            dataAlpha = [];
            a0 = op.a;
            a_current = a0;
            while or(memClose>10*memFar,memClose < memFar/10)
                % Find out if memory usage is ok
                lambdaGuess = (memFar/memClose)^(1/(alpha+2));
                a_current = lambdaGuess*a0;
                Nxi = op.Nxi(a_current);
                NCI = op.NCI(a_current);
                memClose = NCI*16;
                memFar = Nxi*24;
                dataAlpha = [dataAlpha,log(NCI/NCI0)/log(a_current/a0)];%#ok
                alpha = mean(dataAlpha);
            end
            op = op.change_a(a_current);
        end
        function[op] = balanceTime(op)
            Vtest = randn(size(op,2),1);
            tClose = tic;
            q_close_test = op.concretePart*Vtest;%#ok
            tClose = toc(tClose);
            if tClose < 0.01
                tClose = tic;
                for i = 1:10
                    q_close_test = op.concretePart*Vtest;%#ok
                end
                tClose = toc(tClose)/10;
            end
            tFar = tic;
            q_far_test = op.abstractPart(Vtest);%#ok
            tFar = toc(tFar);
            while or(tClose>5*tFar,tClose < tFar/5)
                if tClose > tFar
                    lambda = 0.75;
                else
                    lambda = 1.25;
                    try
                        MemTest = ones(fix(1.5*op.NCI),1); %#ok
                        clear A;
                    catch
                        break
                    end
                end
                if lambda*op.a > 1
                    return
                end
                op = op.change_a(lambda*op.a);
                tClose = tic;
                q_close_test = op.concretePart*Vtest;%#ok
                tClose = toc(tClose);
                if tClose < 0.01
                    tClose = tic;
                    for i = 1:10
                        q_close_test = op.concretePart*Vtest;%#ok
                    end
                    tClose = toc(tClose)/10;
                end
                tFar = tic;
                q_far_test = op.abstractPart(Vtest);%#ok
                tFar = toc(tFar);
            end
        end
        function[op] = computeFullMatrix(op)
            tTotalAssembling = tic;
            op.timesClose.rangesearch = 0;
            
            op.timesClose.NUFFTcloseField = 0;
            
            X1 = op.x(:,1); X2 = op.x(:,2);
            Y1 = op.y(:,1); Y2 = op.y(:,2);
            NX = length(X1);
            NY = length(Y1);
            tComputeInteractions = tic;
            X1_Y1 = repmat(X1,1,NY) - repmat(Y1',NX,1);
            X2_Y2 = repmat(X2,1,NY) - repmat(Y2',NX,1);
            rXY = sqrt(X1_Y1.^2 + X2_Y2.^2);
            rXY(rXY*op.rMax < 1e-15) = 1e-15/op.rMax;
            op.timesClose.computeInteractions = toc(tComputeInteractions);
            op.concretePart = op.kernel.func(rXY);
%             op.concretePart(or(isnan(op.concretePart),...
%                 isinf(op.concretePart))) = 0;
            op.timeTotalAssembling = toc(tTotalAssembling);
            op.full = true;
            op.timesClose.total = toc(tTotalAssembling);
        end
        
        %% Display
        function[] = disp(this)
            printStarLine;
            fprintf('GENERAL INFO \n')
            fprintf('\nOperator of size %s x %s for kernel \n',num2str(size(this,1)),num2str(size(this,2)));
            disp(this.kernel.dilatation(1/this.rMax));
            if this.full
                printStarLine;
                fprtinf('Operator is FULLY assmbled \n');
            else
                rm = this.rMin;
                rM = this.rMax;
                fprintf('Close operator defined for r in [0, %.2f] \n',rm);
                fprintf('Far operator defined for [%.2f, %.2f] \n',rm,rM);
                fprintf('Number of close interactions : %s\n',num2str(this.NCI));
                fprintf('Number of radial quadrature points : %s\n',num2str(length(this.rq.rho) - 1));
                fprintf('Number of 2D quadrature points : %s\n',num2str(this.Nxi));
            end
        end
        function[] = show(this)
            disp(this);
            showTime(this);
            showPerf(this);
            show(this.rq);
        end
        function[] = showTime(this)
            printStarLine;
            tClose = this.timesClose.total;
            tFar = this.timesFar.total;
            tTot = tClose + tFar;
            if this.full
                fprintf('Operator assembled FULLY in %s s\n',tTot);
            else
                fprintf('OFFLINE PERFORMANCE \n')
                fprintf('\n- Operator assembling : %.2f s  :\n',tTot);
                alinea =  '     ';
                
                fprintf([alinea '- Close operator : %.2f s (%.2f %%): \n'],tClose,tClose/tTot*100);
                if this.timesClose.computingTree == 0
                    fprintf([alinea alinea '* KD-tree : 0 %% (precomputed)\n']);
                else
                    fprintf([alinea alinea '* KD-tree : %.2f %% \n'],this.timesClose.computingTree/tClose*100);
                end
                fprintf([alinea alinea '* Range search : %.2f %% \n'],this.timesClose.rangeSearch/tClose*100);
                fprintf([alinea alinea '* Computation of close interactions : %.2f %% \n'],this.timesClose.computeInteractions/tClose*100);
                fprintf([alinea alinea '* Local correction of the far-field (NUFFT) : %.2f %% \n'],this.timesClose.NUFFTcloseField/tClose*100);
                fprintf([alinea alinea '* Assembling close operator : %.2f %% \n'],this.timesClose.assemble/tClose*100);
                
                
                fprintf([alinea '- Far operator : %.2f s (%.2f %%): \n'],tFar,tFar/tTot*100);
                fprintf([alinea alinea '* Radial quadrature : %.2f %% \n'],this.timesFar.radialQuadrature/tFar*100);
                fprintf([alinea alinea '* 2D quadrature : %.2f %% \n'],this.timesFar.quad2D/tFar*100);
            end
            
        end
        function[] = showPerf(this)
            printStarLine;
            alinea =  '     ';
            
            fprintf('ONLINE PERFORMANCE \n')
            % Memory
            fprintf('\n- Memory allocation : %s \n',this.memorySize)
            if this.full
                fprintf('Operator is fully assembled \n')
            else
                sparsityConst = this.NCI/(size(this,1)*size(this,2));
                fprintf([alinea '- Close operator : %.2f %% (%.2f %% sparsity)\n'],this.NCI*2/(this.Nxi*3 + this.NCI*2)*100,100 - sparsityConst*100);
                fprintf([alinea '- Far operator : %.2f %%\n'],this.Nxi/(this.Nxi + this.NCI)*100)
                if (this.NCI*2 > this.Nxi*3)
                    str_for_a = 'reduce';
                else
                    str_for_a = 'increase';
                end
                if or(2*this.NCI > 10*3*this.Nxi,2*this.NCI < 1/10*3*this.Nxi)
                    fprintf('(for more balanced memory allocation, %s parameter ''a'') \n',str_for_a);
                end
            end
            
            % Chrono of MV prod
            if this.full
            else
                if (size(this,1)*size(this,2)<500000)
                    nIt = 20;
                elseif size(this,1)*size(this,2)<5000000
                    nIt = 5;
                else
                    nIt = 1;
                end
                tMV = 0;
                tMVclose = 0;
                tMVfar = 0;
                for i = 1:nIt
                    V = randn(this.N2,1);
                    tMVTic = tic;
                    q = this*V; %#ok
                    tMV = tMV + toc(tMVTic)/nIt;
                    tMVfarTic = tic;
                    q = this.q2d.conv(this.x,this.y,V); %#ok
                    tMVfar = tMVfar + 1/nIt*toc(tMVfarTic);
                    tMVcloseTic = tic;
                    q = this.concretePart*V; %#ok
                    tMVclose = tMVclose + 1/nIt*toc(tMVcloseTic);
                end
            end
            fprintf([alinea '- MV product time : %s s \n'],num2str(tMV));
            fprintf([alinea alinea ' %s %% close, %s %% far \n\n' ],num2str(tMVclose/(tMVfar + tMVclose)*100),num2str(tMVfar/(tMVfar + tMVclose)*100));
        end
        function[str] = memorySize(this)
            NumBytes = this.Nxi*24 + this.NCI*16; % w_nu double, xi_nu(:,1)
            % double, xi_nu(:,2) double, complex numbers in close
            str = Bytes2str(NumBytes);
        end
        %% Others
        function[err,info] = validate(this,V,q,tVal,tol)
            
            if this.full
                disp('Operator is full, nothing to validate');
                return
            end
            if nargin == 1
                V = randn(size(this,2),1);
                V = V/norm(V,1);
                q = this*V;
                tVal = 2;
                tol = this.tol;
            end
            totalTime = tic;
            yy = this.y;
            fun = this.kernel.func;
            
            xx = this.x(1,:);
            rrxy = sqrt((xx(1)-yy(:,1)).^2 +(xx(2)-yy(:,2)).^2);
            tic; ignoredValue = sum(fun(rrxy).*V); t = toc;  %#ok just to count time
            NBEM = min(fix(tVal/t/5),size(this,1));
            if NBEM == 0
                err = NaN;
                warning('Could not validate any sample on such a short amount of time')
            else
                shuffle = randperm(size(this,1));
                q2 = zeros(size(this,1),1);
                beginning = tic;
                try
                    for i = 1:size(this,1)
                        id = shuffle(i);
                        xx = this.x(id,:);
                        rrxy = sqrt((xx(1)-yy(:,1)).^2 +(xx(2)-yy(:,2)).^2);
                        eval = fun(rrxy);
                        eval(or(isinf(eval),isnan(eval))) = 0;
                        q2(i + max(toc(beginning)-tVal,0)*1i) = sum(eval.*V);
                    end
                catch
                    %raised on purpose
                    i = i-1;
                end
                valSample = shuffle(1:i);
                err = norm(q(valSample)-q2(1:i),'inf');
                t = toc(totalTime);
                if err > tol
                    warning('Found samples with greater error than required accuracy')
                    info.badSamples = valSample(abs(q(valSample)-q2(1:i))>tol);
                    info.maxErr = err;
                else
                    propVal = min(i/size(this,1)*100,100);
                    fprintf('\n Validated %s%% of the samples in %s seconds \n',num2str(propVal),num2str(t));
                    info = [];
                    disp(':-)');
                end
                
                
            end
        end
        function[rm] = rMin(this)
            rm = this.a*this.rMax;
        end
        
    end
    methods (Static)
        function[rM] = rMaxCalc(X,Y)
            dist = @(x,x0) max( sqrt( (x(:,1)-x0(1)).^2 + (x(:,2)-x0(2)).^2) );
            rX   = dist(X,mean(X,1));
            rY   = dist(Y,mean(Y,1));
            rXY  = dist(mean(X,1),mean(Y,1));
            rM = rXY + rX  + rY;
        end
        function[aa] = default_a(N1,N2)
            aa = 1/(4*sqrt(sqrt(N1*N2)));
        end
    end
end
classdef (Abstract) singularKernel
    % Kernels for which we can compute explicitly integrals along segments.
    
    properties
        id
        k@Kernel
    end
    methods
        function[s,W] = singularKernelQuadrature(this,X,A,B,N)
            % Given two points A and B in R², a list of points X1,...,XM in R², and a
            % singular kernel K, computes a list of order N Gaussian quadrature rules on
            % the segment [A,B], that is M lists of N points y_1^k, y_2^k,...y_N^k
            % lying in [A,B] and M vectors of weights w_1^k, ..., w_N^k such that
            % for all k, and any function f sufficiently smooth
            % \int_[A,B] K(x^k,y)f(y)dy approx sum_{i} w_i^k f(y_i^k)
            % applying the Golub-Welsch algorithm
            % The problem is rescaled to computing weights and nodes for
            % \int_{0}^{1} K(x^k,A + (B - A)*u)f(A + (B -A)*u)du
            % = sum_{i} wtilde_i^k f(A + (B - A) utilde_i^k)
            % The first integral is recovered by
            % int_{[A,B]}K(x^k,y)f(y)dy = norm(B-A) \int_{0}^1 K(x^k,A +
            % (B-A)u) f(A + (B-A)u) du  approx norm(B-A) \sum_{i}wtilde_i^k f(A + (B-A)utilde_i^k)
            % Thus
            % s_i^k = utilde_i^k %  It is precisely the abscissa on [A,B]
            % w_i^k = norm(B-A)*wtilde_i^k
            
            if N > Nmax(this)
                warning(['It is not possible to compute a %s-order quadrature,' ...
                    'using Nmax = %s instead'],num2str,(N),num2str(this.Nmax))
                N = Nmax(this);
            end
            M = size(X,1);
            l = norm(B-A);
            
            % Prepare outputs
            W = zeros(N,M);
            s = zeros(N,M);
            
            a = zeros(N,M);
            b = zeros(N,M);
            ints = zeros(2*N,M);
            orthoPols = cell(N,1);
            orthoPols{1} = ones(1,M);
            ints(1,:) = this.IN(X,A,B,0).';
            ints(2,:) = this.IN(X,A,B,1).';
            prpr = ints(1,:);
            mu0  = prpr;
            xprpr = ints(2,:);
            a(1,:) = xprpr./prpr;
            b(1,:) = 0;
            for r = 2:N
                if r==2
                    orthoPols{r} = [zeros(1,M);orthoPols{r-1}] - ...
                        [repmat(a(r-1,:),r-1,1).*orthoPols{r-1};zeros(1,M)];
                else
                    orthoPols{r} = [zeros(1,M);orthoPols{r-1}] - ...
                        [repmat(a(r-1,:),r-1,1).*orthoPols{r-1};zeros(1,M)]-...
                        [repmat(b(r-1,:),r-2,1).*orthoPols{r-2};zeros(2,M)];
                end
                pr_1pr_1 = prpr;
                ints(2*r-1,:) = this.IN(X,A,B,2*r-2).';
                ints(2*r,:) = this.IN(X,A,B,2*r-1).';
                sqP = squarePol(orthoPols{r});
                Psquare = [sqP;zeros(1,M)];
                xPsquare = [zeros(1,M);sqP];
                prpr = sum(Psquare.*ints(1:2*r,:),1);
                xprpr = sum(xPsquare.*ints(1:2*r,:),1);
                a(r,:) = xprpr./prpr;
                b(r,:) = prpr./pr_1pr_1;
            end
            for m = 1:M
                am = a(:,m);
                bm = b(:,m);
                mu0m = mu0(m);
                Jm = diag(sqrt(bm(2:end)),-1) + diag(am) + diag(sqrt(bm(2:end)),1);
                [Pm,Dm] = eig(Jm);
                sm = diag(Dm);
                omegam = mu0m*Pm(1,:).^2;
                omegam = omegam(:);
                s(:,m) = sm;
                W(:,m) = l*omegam;
            end
        end
        
        function[I] = I0seg(this,X,A,B)
            % returns \int_{[A,B]} k(X,Y)dY
            I = this.IN(X,A,B,0);
        end
        
    end
    methods (Static)
        function [a,b,d] = parameters_singInt(X,A,B)
            
            N = normalVector(A,B);
            AB = B-A;
            XB = B - X;
            XA = A - X;
            l = sqrt(cWise_dot(AB,AB));
            u = AB./[l l];
            a = cWise_dot(XA,u);
            b = cWise_dot(XB,u);
            d = abs(cWise_dot(N,XA));
            
            
            
            function[Z] = cWise_dot(X,Y)
                Z = X(:,1).*Y(:,1) + X(:,2).*Y(:,2);
            end
            
        end
    end
    
    methods (Abstract)
        I = IN(this,x,A,B,N)
        % I = IN(this,x,A,B,N) :
        % computes the integral I = \int_{0}^1 k(X,A + (B-A)*u) u^n du
        Nm = Nmax(this)
%         I = Icomp_r(this,r,t0,rt0_ders,A,B)
        % computes \int_{a}^b k(r(t0), r(t)) dt
        F = primitiveOfKernel(this,s)
    end
    
end


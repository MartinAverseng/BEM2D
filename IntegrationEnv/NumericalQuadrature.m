classdef NumericalQuadrature
    
    properties
        num@double
        order@double
        xhat@double
        what@double        
    end
    
    methods
        function[this] = NumericalQuadrature(n)            
            this.num = n;
            [x,w,o] = NumericalQuadrature.quadRule(this.num);
            this.order = o;
            this.xhat = x;
            this.what = w;            
        end
        function[] = disp(this)
            fprintf('%d-point Quadrature rule of order %d on the segment [0,1]\n\n',this.num,this.order)
            x = this.xhat;
            w = this.what;
            T = table(x,w);
            disp(T);
        end
        function[] = show(this)
            figure;            
            plot(this.xhat,this.what,'*');
            title(sprintf('Quadrature rule of order %d on [0,1]',this.order));
            xlim([0,1]);
            set(gca,'XTick',[0;this.xhat;1]);            
            set(gca,'YTick',unique(this.what));
            xlabel('$\hat{x}$','Interpreter','latex','FontSize',16)
            ylabel('$\hat{w}$','Interpreter','latex','FontSize',16)
            set(gca,'FontSize',12);
            
        end
        function[N] = nint(this,mesh)
            N = length(this.weights(mesh));
        end
        function[] = showOn(this,mesh)
            mesh.plot;
            hold on
            x = this.intPoints(mesh);            
            plot(x(:,1),x(:,2),'*');
        end
        function[sAB,XAB,wAB] = transport(this,sA,sB,A,B)
            sAB = sA + this.xhat*(sB - sA);
            XAB1 = A(1) + this.xhat*(B(1) - A(1));
            XAB2 = A(2) + this.xhat*(B(2) - A(2));
            XAB = [XAB1 XAB2];
            wAB = this.what*(sB-sA);
        end
        function[x] = intPoints(this,mesh)
            n = this.num;
            x = zeros(mesh.nseg*n,2); 
            [zx1,zx2,zy1,zy2] = mesh.edgesCoords;            
            for i = 1:n
                x(i:n:end,1) = zx1 + (zx2-zx1)*this.xhat(i);
                x(i:n:end,2) = zy1 + (zy2-zy1)*this.xhat(i);                
            end
        end
        function[w] = weights(this,mesh)
            n = this.num;
            w = zeros(mesh.nseg*n,1);
            X = mesh.edgesCoords;
            zx1 = X(:,1); zx2 = X(:,3); zy1 = X(:,2); zy2 = X(:,4);
            for i = 1:n               
                w(i:n:end) = sqrt((zx1-zx2).^2 + (zy1 - zy2).^2)*this.what(i);
            end
        end
    end
    
    methods(Static)
        function[x,w,o] = quadRule(n)
           o = 2*n - 1;
           [x,w] = lgwt(n,0,1);
        end
    end
end


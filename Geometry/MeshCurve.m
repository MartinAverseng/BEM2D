classdef MeshCurve < handle
    properties
        vertices; %
        segments; % The list of couples (i,j) such that [t(i),t(j)] is a
        %segment of the mesh, where t is the list of vertices
        curve@SimpleCurve;
        tVertices;
        repartition = []; %optional
        bounds = [-1,1];
        edgesCoords;
        sVertices;
        length
        L
    end
    
    methods
        function[this] = MeshCurve(c,N,rep,bds)
            % N is the number of segments.
            this.curve = c;
            l = length(c);
            s = linspace(0,l,N+1); 
            s = s(2:end-1);
            a = c.I(1);
            b = c.I(2);
            if ~ exist('rep','var')||isempty(rep)
                this.repartition = [];
                t = c.s_1(s);
                t = [a,t,b]';
                this.bounds = [c.I(1),c.I(2)];
            else
                this.repartition = rep;
                this.bounds = bds;
                alpha = bds(1);
                beta = bds(2);
                assert(abs(rep(alpha)-a)<1e-10);
                assert(abs(rep(beta)-b)<1e-10);
                t = rep(linspace(alpha,beta,N+1)');
            end
            
            % When the curve is closed, we have M(t1) = M(tN)
            if c.closed
                assert(N>1);% otherwise, the curve cannot be closed !
                t = t(1:end-1);
                N = N-1;
            end
            s = (t - a)/(b-a)*l;
            t = unique(c.s_1(s));
            this.tVertices = t(:);
            N = length(t)-1;
            x = c.x(t);
            y = c.y(t);
            this.vertices = [x,y];
            this.segments = [(1:N)',(2:N+1)'];
            if c.closed
                this.segments = [this.segments;N+1,1];
            end
            xA = this.vertices(this.segments(:,1),1);
            yA = this.vertices(this.segments(:,1),2);
            xB = this.vertices(this.segments(:,2),1);
            yB = this.vertices(this.segments(:,2),2);
            this.edgesCoords = [xA,yA,xB,yB];
            this.length = sqrt((xB-xA).^2 + (yB-yA).^2);
            this.L = sum(this.length);
            this.sVertices = [0;cumsum(this.length)];
        end
        function[this] = remesh(this,N)
            this = MeshCurve(this.curve,N,this.repartition,this.bounds);
        end
        function[N] = nvert(this)
            N = size(this.vertices,1);
        end
        function[N] = nseg(this)
            N = size(this.segments,1);
        end
        function[] = plot(this,withVecs)
            if nargin==1
                withVecs = true;
            end
            X = this.edgesCoords;
            x1 = X(:,1); y1 = X(:,2); x2 = X(:,3); y2 = X(:,4);
            
            X = [x1(:), x2(:)]';
            Y = [y1(:), y2(:)]';
            X = X(:);
            Y = Y(:);
            plot(X,Y);
            axis equal;
            
            hold on
            mid = this.middle;
            plot(mid(:,1),mid(:,2),'o','MarkerSize',3);
            if withVecs
                t = this.tangentVector;
                n = this.normalVector;
                quiver(mid(:,1),mid(:,2),n(:,1),n(:,2),0);
                quiver(mid(:,1),mid(:,2),t(:,1),t(:,2),0);
            end
        end
        function[vert] = vertex(this,i)
            verts = this.vertices;
            vert = verts(i,:);
        end
        function[seg] = segment(this,i)
            segs = this.segments;
            seg = segs(i,:);
        end
        function[tl] = Tl_t(this)
            % returns the list of applications Tl : [0,1] -> [tl,tl+1]
            % defined by Tl(u) = a*u + b
            a = diff(this.tVertices);
            b = this.tVertices(1:end-1);
            tl = [a(:),b(:)];
            
        end
        function[MM] = M(this,s)
            s = s(:);
            X = this.edgesCoords;
            x1 = X(:,1); y1 = X(:,2); x2 = X(:,3); y2 = X(:,4);
            ids = this.whichSegment_s(s);
            sVert = this.sVertices;
            sloc = (s - sVert(ids))./(sVert(ids+1)-sVert(ids));
            x = x1(ids,:) + sloc.*(x2(ids)-x1(ids));
            y = y1(ids,:) + sloc.*(y2(ids)-y1(ids));
            MM = [x(:),y(:)];
        end
        function[ss] = s(this,t)
            % Returns the curvilinear abscissa
            l = this.length;
            int_l = [0;cumsum(l(:))];
            ids = this.whichSegment_t(t);
            tVert = this.tVertices;
            tloc = (t - tVert(ids))./(tVert(ids+1)-tVert(ids));
            ss = int_l(ids) + tloc.*l(ids);
        end
        function[sA,sB] = s_seg(this,seg_num)
            sVert = this.sVertices;
            sA = sVert(seg_num);
            sB = sVert(seg_num + 1);
        end
        function[ids] = whichSegment_s(this,s)
            assert(sum((s-this.s(this.curve.I(2)))>1e-11)==0)
            assert(sum(s<-1e-11)==0)
            s = s(:);
            sVert = this.sVertices;
            ids = sum(repmat(s,1,length(sVert)-1) - repmat(sVert(1:end-1)',length(s),1) >= 0,2);
        end
        function[ids] = whichSegment_t(this,t)
            assert(sum((t-this.curve.I(2))>1e-11)==0)
            assert(sum((t-this.curve.I(1))<-1e-11)==0)
            t = t(:);
            tVert = this.tVertices;
            ids = sum(repmat(t,1,length(tVert)-1) - repmat(tVert(1:end-1)',length(t),1) >= 0,2);
        end
        function[mid] = middle(this,i)
            if nargin ==1
                i = 1:this.nseg;
            end
            sVert = this.sVertices;
            ds = diff(sVert);
            mid = this.M(sVert(i) + ds(i)/2);
        end
        
        function[t] = tangentVector(this,i)
            if nargin==1
                i = 1:this.nseg;
            end
            X = this.edgesCoords;
            x1 = X(:,1); y1 = X(:,2); x2 = X(:,3); y2 = X(:,4);
            t = [x2(i)-x1(i),y2(i)-y1(i)];
        end
        function[n] = normalVector(this,i)
            if nargin==1
                i = 1:this.nseg;
            end
            t = this.tangentVector(i);
            directRot = [0 -1; 1 0];
            switch this.curve.boundedSide
                case 'right'
                    Rot = directRot;
                case 'left'
                    Rot = directRot';
                case 'none'
                    Rot = directRot;
            end
            n = (Rot*t')';
        end
        function[] = animate(this)
            figure
            this.plot;
            hold on;
            N = 200;
            t = linspace(this.curve.I(1),this.curve.I(2),N);
            Mt = this.M(t);
            h2 = plot(Mt(1,1),Mt(1,2),'k*');
            for i = 1:N
                Mi = Mt(i,:);
                delete(h2);
                h2 = plot(Mi(1),Mi(2),'k*');
                pause(0.01);
                drawnow;
            end
            hold off
        end        
        function[this] = refine(this,localLength)
            warning('not supported yet');
        end
    end
    
end


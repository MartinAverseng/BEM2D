classdef SimpleCurve
    % Parametric curve, can be closed but must not intersect itself
    
    properties
        x; % Function of 1 argument (the parameter t)
        y; %
        I; % interval of the values for the parameter t
        closed@logical;
        boundedSide;
    end
    
    methods
        function[this] = SimpleCurve(xx,yy,II,bS)
            this.x = xx;
            this.y = yy;
            this.I = II;
            a = II(1); b = II(2); Ma = [xx(a);yy(a)]; Mb = [xx(b);yy(b)];
            this.closed = norm(Ma-Mb)<1e-11;
            if this.closed
                assert(logical(exist('bS','var')),'this curve is closed. You must pass the boundedSide (value "left" or "right") in argument');
                assert(ismember(bS,{'left','right'}),['This curve is closed. the value used for argument bS is incorrect. Please use one of the choices : "left" or "right". \n' ...
                     'choose left if the bounded component of the plane lies at the left of the curve, and right otherwise.'])
                this.boundedSide = bS;
            else
                this.boundedSide = 'none';
                
            end
        end
        function[] = plot(this,varargin)
            N = 100;
            II = this.I;
            t = linspace(II(1),II(2),N);
            xt = this.x(t);
            yt = this.y(t);
            plot(xt,yt,varargin{:});
            axis equal
        end
        function[S] = s(this,t)
            Integrand = @(t)(sqrt(this.dx(t).^2 + this.dy(t).^2));
            [~,i] = sort(t);
            S = zeros(size(t));
            S(i(1)) = integral(Integrand,this.I(1),t(1));
            for l = 2:length(t)
                S(i(l)) = S(i(l-1)) + integral(Integrand,t(i(l-1)),t(i(l)));
            end
        end
        function[DX] = dx(this,t)
            DX = 1e6*(this.x(t + 1e-6) - this.x(t));
        end
        function[DY] = dy(this,t)
            DY =  1e6*(this.y(t + 1e-6) - this.y(t));
        end
        function[r] = length(this)
            dx = @(t)(1e6*(this.x(t + 1e-6) - this.x(t)));
            dy = @(t)(1e6*(this.y(t + 1e-6) - this.y(t)));
            r = integral(@(t)(sqrt(dx(t).^2 + dy(t).^2)),this.I(1),this.I(2));
        end
        function[t] = s_1(this,s)
            [~,Idx] = sort(s);
            t = zeros(size(s));
            t0 = -1;
            for i = 1:length(Idx)
                t(Idx(i)) = fzero(@(t)(this.s(t) - s(Idx(i))),t0);
                t0 = t(Idx(i));
            end
        end
    end
end


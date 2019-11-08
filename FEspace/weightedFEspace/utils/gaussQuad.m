function [ x,omega ] = gaussQuad( w,N,alpha,beta )
% finds a N points Gaussian quadrature for function w using Golub-Welsch
% algorithm.

%wprime = @(u)(w(alpha + (beta - alpha)*u));
wprime = w;
a = zeros(N,1);
b = zeros(N,1);
p = cell(N,1);
p{1} = @(x)1;
prpr = integral(@(x)(wprime(x).*p{1}(x).^2),alpha,beta);
mu0 = prpr;
xprpr = integral(@(x)(x.*wprime(x).*p{1}(x).^2),alpha,beta);
a(1) = xprpr/prpr;
b(1) = 0;

for r = 2:N
    if r==2
        p{r} = @(x)((x-a(r-1)).*p{r-1}(x));
    else
        p{r} = @(x)((x-a(r-1)).*p{r-1}(x) - b(r-1)*p{r-2}(x));
    end
    pr_1pr_1 = prpr;
    prpr = integral(@(x)((p{r}(x)).^2.*wprime(x)),alpha,beta);
    xprpr = integral(@(x)(x.*(p{r}(x)).^2.*wprime(x)),alpha,beta);
    a(r) = xprpr/prpr;
    b(r) = prpr/pr_1pr_1;
end

J = diag(sqrt(b(2:end)),-1) + diag(a) + diag(sqrt(b(2:end)),1);
[P,D] = eig(J);
s = diag(D);
%x = alpha + (beta - alpha)*s;
x = s;
omega = mu0*P(1,:).^2;
omega = omega(:);
%omega = (beta-alpha)*omega(:);

end


clear
Nx = 36; Nk = 20;
xs = 1; xf = 2*pi;
xx = linspace(xs,xf,Nx);
x = 2*pi*(xx-xs)/(xf-xs) - pi;

kk = (-Nk/2:Nk/2);

c = fourierC(Nk);
ff = zeros(1,Nx);

for jj=1:Nx    
   ff(jj) = sum(c .*  exp(1i*kk*x(jj)) );
end


plot(xx,ff,xx,f(x),'.')

% a0 = 1/2/pi .* integral(f,-pi,pi);
% af = @(n) 1/pi .* integral(@(x)cos(n.*x).*f(x),-pi,pi);
% bf = @(n) 1/pi .* integral(@(x)sin(n.*x).*f(x),-pi,pi);

% a = zeros(1,Nk);
% b = zeros(1,Nk);
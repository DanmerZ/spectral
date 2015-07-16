clear
N = 6;
syms x f;
f = @(x) 1./(1+x.^2);
X = linspace(-5,5,N);
F = f(X);

xx = linspace(-5,5,100);
ff = f(xx);

C = [];
for ii=1:N
    C = [C; @(x) ]
end


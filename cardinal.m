clear;

N = 50;
NN = 100;

xmin = -2;
xmax = 2;

% x = linspace(xmin,xmax,N);
% xx = linspace(xmin,xmax,NN);
x = linspace(0,2*pi*(N-1)/N,N);
xx = linspace(0,2*pi*(N-1)/N,NN);

h = x(2) - x(1);

y = 1.0 ./ (1.0 + x.*x);
yy = xx;

for k=1:NN
    yy(k) = 0;
    for j=1:N
%         sinc = sin(pi*(xx(k)-(j-N/2)*h)/h) / (pi*(xx(k)-(j-N/2)*h)/h);
        sinc = sin(N * (xx(k) - x(j))/2) * cot((xx(k) - x(j))/2) / N;
        yy(k) = yy(k) + y(j)*sinc;
    end
end

plot(x,y,'o',xx,yy,'-')

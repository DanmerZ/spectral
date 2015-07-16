clear;

N = 40;
NN = 100;

xmin = -1;
xmax = 1;

x = linspace(xmin,xmax,N);
for i=1:N
     x(i) = -cos( (2*i-1)*pi / (2*(N+1)) );
end
xx = linspace(xmin,xmax,NN);

y = 1 ./ (1 + x.*x);
yy = xx;

for k=1:NN
    yy(k) = 0;
    for i=1:N
        c = 1;
        for j=1:N
            if i ~= j
                c = c * (xx(k) - x(j)) / (x(i) - x(j));
            end
        end
        yy(k) = yy(k) + y(i) * c; 
    end    
end

% plot(x,y)
plot(x,y,'ko',xx,yy,'r')
clear
N = 240;
x = linspace(-pi,pi,N)';
h = x(2) - x(1);

% fx = -1.5*cos(x) - 4.5*cos(2*x);
fx = exp(-1 + cos(x));

nn = 4;
phi = [];
for jj=0:nn-1
    phi = [phi cos(jj*x)];
end

H = zeros(nn);
f_ = [];
for i = 1:nn
   for j=1:nn
       H(i,j) = trapz(x,phi(:,i).*H_(phi(:,j),h,x));
   end
   f_ = [f_; trapz(x,fx.*phi(:,i))];
end


% f_ = [trapz(x,fx.*cos(x)); trapz(x,fx.*cos(2*x))];

a = H\f_;

u = phi * a;
plot(x,u)


function [ f ] = f(x)
% f = sin(x).* (x+pi < pi*ones(1,length(x)));
% f = x
% f = exp(-cos(x).^2/sin(x).^2);
p = 0.1;
f = (1-p^2)./(1+p^2-2*p*cos(x));
end

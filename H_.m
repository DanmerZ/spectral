function [ H_ ] = H_( f,h,x )
%H_ computes internal H phi_i

d1 = gradient(f,h);
d2 = gradient(d1,h);

% H_ = d2 - .5*f;
H_ = d2 + (cos(x) + cos(x).^2).*f;


end


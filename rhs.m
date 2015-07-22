function [res] = rhs( uf,alpha,ikx2,iky2 )
%rhs Right hand side of equation alpha*(uxx+uyy) or etc
[uxx,uyy] = sd2_2d(uf,ikx2,iky2);
res = alpha*(uxx+uyy);


end


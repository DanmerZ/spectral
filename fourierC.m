function [ c ] = fourierC( Nk )

cf = @(n) .5/pi * integral(@(x) exp(-1i*n.*x).*f(x),-pi,pi);

kk = (-Nk/2:Nk/2);
c = zeros(1,Nk);

for jj=kk
    c(jj+Nk/2+1) = cf(jj);    
end

end


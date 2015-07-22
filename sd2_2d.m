function [ uxx, uyy ] = sd2_2d( uf, ikx2, iky2 )
%sd2_2d Second spectral derivatives in 2D

uxx = real(ifft2(ikx2.*uf));
uyy = real(ifft2(iky2.*uf));


end


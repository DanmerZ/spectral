function [ ux,uy ] = sd1_2d(uf,ikx,iky)
%sd1_2d first spectral derivatives

ux = real(ifft2(ikx.*uf));
uy = real(ifft2(iky.*uf));


end


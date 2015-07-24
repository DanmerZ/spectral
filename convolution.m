clear all;
clc;
close all;

NX = 512;
dx = 2*pi/NX;

i = meshgrid(1:NX,1);
x = i*dx;

kx = [0:NX/2 -NX/2+1:-1];

f = exp(-(x-pi).^2);
ff = f.*f;
dffx = -4*(x-pi).*ff;

ffhat = fft(ff);

ffhat2 = ifft( fft(fft(f)).^2 );



plot(x,-dffx+real(ifft(1i*kx.*ffhat2))/NX)


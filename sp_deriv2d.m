clear all;

% resolution in x and y 
NX = 128; NY = 128;
dx = 2*pi/NX; dy = 2*pi/NY;

[i,j] = meshgrid(1:NX,1:NY);
x = i*dx; y = j*dy;

% matrices of wavenumbers in x and y
kx = 1i*ones(1,NY)'*[0:NX/2 -NX/2+1:-1];            %ones(1,NY)'*fftshift(-NX/2+1:NX/2);
ky = 1i*[0:NY/2 -NY/2+1:-1]'*ones(1,NX);            %fftshift([-NY/2+1:NY/2]')*ones(1,NX);

f = exp(-((i*dx-pi).^2 + (j*dy-pi).^2));
% f = sin(-((x-pi).^2 + (y-pi).^2));

% 1st analytical derivatives
fx = -2*(x-pi).*f;
fy = -2*(y-pi).*f;
% fx = -2*(x-pi).*cos(-((x-pi).^2 + (y-pi).^2));
% fy = -2*(y-pi).*cos(-((x-pi).^2 + (y-pi).^2));

% 2nd analytical derivatives
fxx = -2*f.*(1-2*(x-pi).^2);
fyy = -2*f.*(1-2*(y-pi).^2);
% fxx = -2*cos(-((x-pi).^2 + (y-pi).^2)) + 4*f.*(x-pi).^2;
% fyy = -2*cos(-((x-pi).^2 + (y-pi).^2)) + 4*f.*(y-pi).^2;

% Mixed derivative
fxy = 4*(x-pi).*(y-pi).*f;

% 2D Fourier transform
f_hat = fft2(f);

% 1st spectral derivatives - OK           f' -> ik fft(f) 
fs_x = real(ifft2(kx.*f_hat));
fs_y = real(ifft2(ky.*f_hat));

% 2nd spectral derivatives - NOT OK       f'' -> (ik)^2 fft(f)
fs_xx = real(ifft2(kx.*kx.*f_hat));
fs_yy = real(ifft2(ky.*ky.*f_hat));

% spectral mixed derivative
fs_xy = real(ifft2(kx.*ky.*f_hat));

subplot(2,2,1); mesh(fs_xx); title('Spectral 2nd deriv.');
subplot(2,2,3); mesh(fxx); title('Analytical 2nd deriv.');
subplot(2,2,2); mesh(fs_xy); title('Spectral mixed deriv.');
subplot(2,2,4); mesh(fxy); title('Analytical mixed deriv.');

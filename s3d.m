clear all;

NX = 32; NY = 32; NZ = 32;
dx = 2*pi/NX; dy = 2*pi/NY; dz = 2*pi/NZ; 


[i,j,k] = meshgrid(1:NX,1:NY,1:NZ);
x = i*dx; y = j*dy; z = k*dz;

f = exp(-((x-pi).^2 + (y-pi).^2 + (z-pi).^2));
% f = sin(-((x-pi).^2 + (y-pi).^2 + (z-pi).^2));

fx = -2*(x-pi).*f;
fy = -2*(y-pi).*f;
fz = -2*(z-pi).*f;

% fx = -2*(x-pi).*cos(-((x-pi).^2 + (y-pi).^2 + (z-pi).^2));
% fy = -2*(y-pi).*cos(-((x-pi).^2 + (y-pi).^2 + (z-pi).^2));
% fz = -2*(z-pi).*cos(-((x-pi).^2 + (y-pi).^2 + (z-pi).^2));

fxx = -2*f.*(1-2*(x-pi).^2);
fyy = -2*f.*(1-2*(y-pi).^2);
fzz = -2*f.*(1-2*(z-pi).^2);
% 
% fxx = -2*cos(-((x-pi).^2 + (y-pi).^2 + (z-pi).^2)) + 4*f.*(x-pi).^2;
% fyy = -2*cos(-((x-pi).^2 + (y-pi).^2 + (z-pi).^2)) + 4*f.*(y-pi).^2;
% fzz = -2*cos(-((x-pi).^2 + (y-pi).^2 + (z-pi).^2)) + 4*f.*(z-pi).^2;


kx = [0:NX/2 -NX/2+1:-1];
ky = [0:NY/2 -NY/2+1:-1];
kz = [0:NZ/2 -NZ/2+1:-1];

ikx = zeros(NX,NY,NZ);
iky = zeros(NX,NY,NZ);
ikz = zeros(NX,NY,NZ);

for ii=1:NZ
    ikx(:,:,ii) = 1i*ones(1,NY)'*kx;
    iky(:,:,ii) = 1i*ky'*ones(1,NX);
end
for ii=1:NX
    ikz(ii,:,:) = 1i*ones(1,NY)'*kz;
end

f_hat = fftn(f);

fs_x = real(ifftn(ikx.*f_hat));
fs_y = real(ifftn(iky.*f_hat));
fs_z = real(ifftn(ikz.*f_hat));


fs_xx = real(ifftn(ikx.*ikx.*f_hat));
fs_yy = real(ifftn(iky.*iky.*f_hat));
fs_zz = real(ifftn(ikz.*ikz.*f_hat));

% 
 subplot(1,2,1); scatter3(x(:),y(:),z(:),15,f(:),'fill'), colorbar
 subplot(1,2,2); slice(x,y,z,fxx-fs_xx,pi,pi,pi), colorbar
% subplot(1,3,1); slice(x,y,z,imag(ikx),pi,pi,pi)
% subplot(1,3,2); slice(x,y,z,imag(iky),pi,pi,pi)
% subplot(1,3,3); slice(x,y,z,imag(ikz),pi,pi,pi)

% subplot(1,3,1); slice(x,y,z,fs_x,pi,pi,pi)
% subplot(1,3,2); slice(x,y,z,fs_y,pi,pi,pi)
% subplot(1,3,3); slice(x,y,z,fs_z,pi,pi,pi)

% colormap('Jet'), colorbar
% subplot(1,3,1); slice(x,y,z,fs_xx,pi,pi,pi)
% subplot(1,3,2); slice(x,y,z,fs_yy,pi,pi,pi)
% subplot(1,3,3); slice(x,y,z,fs_zz,pi,pi,pi)
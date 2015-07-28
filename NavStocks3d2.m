clear all;
clc;
close all;

nu = 1.0; % viscosity

NX = 32; NY = 32; NZ = 32;
dx = 2*pi/NX; dy = 2*pi/NY; dz = 2*pi/NZ;

[x,y,z] = meshgrid(dx*[1:NX],dy*[1:NY],dz*[1:NZ]);

% components of velocity
% u = U_x, v = U_y, w = U_z

u = exp(-((x-pi).^2+(y-pi).^2+(z-pi).^2));
v = zeros(NX,NY,NZ);
w = zeros(NX,NY,NZ);

% SPECTRAL
% wave numbers
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

ikx2 = ikx.*ikx;
iky2 = iky.*iky;
ikz2 = ikz.*ikz;

k2 = abs(ikx.^2 + iky.^2 + ikz.^2);

% Fourier transform of velocity components

uf = fftn(u); uff = fftn(uf);
vf = fftn(v); vff = fftn(vf);
wf = fftn(w); wff = fftn(wf);

% Fourier transform of advection term
% as convolution

advxf = ikx.*uff.*uff + iky.*vff.*uff + ikz.*wff.*uff;
advyf = ikx.*uff.*vff + iky.*vff.*vff + ikz.*wff.*vff;
advzf = ikx.*uff.*wff + iky.*vff.*wff + ikz.*wff.*wff;

advx = real(ifftn(advxf))/NX;

% subplot(2,2,1); slice(x,y,z,real(ifftn(advxf)),pi,pi,pi),colorbar
% subplot(2,2,2); slice(x,y,z,real(ifftn(advyf)),pi,pi,pi),colorbar
% subplot(2,2,3); slice(x,y,z,real(ifftn(advzf)),pi,pi,pi),colorbar
% subplot(2,2,4); slice(x,y,z,u,pi,pi,pi),colorbar

surf(advx(:,:,NZ))


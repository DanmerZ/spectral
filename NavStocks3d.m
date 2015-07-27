clear all;
clc;
close all;

nu = 1.0; % viscosity

NX = 16; NY = 16; NZ = 16;
dx = 2*pi/NX; dy = 2*pi/NY; dz = 2*pi/NZ;

[x,y,z] = meshgrid(dx*[1:NX],dy*[1:NY],dz*[1:NZ]);

% components of velocity
% u = U_x, v = U_y, w = U_z

u = exp(-((x-pi).^2+(y-pi).^2+(z-pi).^2));
v = zeros(NX,NY,NZ);
w = zeros(NX,NY,NZ);

% analytical derivatives

ux = -2.0*(x-pi).*u; uy = -2.0*(y-pi).*u; uz = -2.0*(z-pi).*u;
vx = 0.0; vy = 0.0; vz = 0.0;
wx = 0.0; wy = 0.0; wz = 0.0;

uxx = -2.0*u.*(1-2.0*(x-pi).^2);
uyy = -2.0*u.*(1-2.0*(y-pi).^2);
uzz = -2.0*u.*(1-2.0*(z-pi).^2);

vxx = 0.0;
vyy = 0.0;
vzz = 0.0;

wxx = 0.0;
wyy = 0.0;
wzz = 0.0;

% analytical advection term derivatives

advx = u.*ux + v.*uy + w.*uz;
advy = u.*vx + v.*vy + w.*vz;
advz = u.*wx + v.*wy + w.*wz;

% analitycal dissipation term

dissx = nu*(uxx + uyy + uzz);
dissy = nu*(vxx + vyy + vzz);
dissz = nu*(wxx + wyy + wzz);

% analytical Poisson pressure resolving
% p(1,1,1) = 0

p = zeros(NX,NY,NZ);



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

uf = fftn(u);
vf = fftn(v);
wf = fftn(w);

% Spectral velocity derivatives

us_x = real(ifftn(ikx.*uf));
us_y = real(ifftn(iky.*uf));
us_z = real(ifftn(ikz.*uf));

us_xx = real(ifftn(ikx2.*uf));
us_yy = real(ifftn(iky2.*uf));
us_zz = real(ifftn(ikz2.*uf));

vs_x = real(ifftn(ikx.*vf));
vs_y = real(ifftn(iky.*vf));
vs_z = real(ifftn(ikz.*vf));

vs_xx = real(ifftn(ikx2.*vf));
vs_yy = real(ifftn(iky2.*vf));
vs_zz = real(ifftn(ikz2.*vf));

ws_x = real(ifftn(ikx.*wf));
ws_y = real(ifftn(iky.*wf));
ws_z = real(ifftn(ikz.*wf));

ws_xx = real(ifftn(ikx2.*wf));
ws_yy = real(ifftn(iky2.*wf));
ws_zz = real(ifftn(ikz2.*wf));

% Spectral advection term

advs_x = u.*us_x + v.*us_y + w.*us_z;
advs_y = u.*vs_x + v.*vs_y + w.*vs_z;
advs_z = u.*ws_x + v.*ws_y + w.*ws_z;

% Spectral dissipative term

diss_x = nu*(us_xx + us_yy + us_zz);
diss_y = nu*(vs_xx + vs_yy + vs_zz);
diss_z = nu*(ws_xx + ws_yy + ws_zz);

% Spectral pressure Poisson resolving

p_hat = (-ikx.*fftn(advx)-iky.*fftn(advy)-ikz.*fftn(advz))./k2;
p_hat(1,1,1) = 0.0;
p_s = real(ifftn(p_hat));
% prhs = real(ifftn(-ikx.*fftn(advx)))+real(ifftn(-iky.*fftn(advy)))+real(ifftn(-ikz.*fftn(advz)));



% scatter3(x(:),y(:),z(:),15,advx(:)),colorbar
% slice(x,y,z,advx,pi,pi,pi),colorbar
surf(p_s(:,:,NZ/2))

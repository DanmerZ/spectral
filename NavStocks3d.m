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

% finite difference Poisson pressure resolving
% p(1,1,1) = 0

p = zeros(NX,NY,NZ);
p1 = zeros(NX,NY,NZ);
pxx = p1; pyy = p1; pzz = p1;
advx_x = p1; advy_y = p1; advz_z = p1;

for ii=2:NX-1
    for jj=2:NX-1
        for kk=2:NY-1
            pxx(ii,jj,kk)=(p(ii+1,jj,kk)-2*p(ii,jj,kk)+p(ii-1,jj,kk))/(dx.^2);
            pyy(ii,jj,kk)=(p(ii,jj+1,kk)-2*p(ii,jj,kk)+p(ii,jj-1,kk))/(dy.^2);
            pzz(ii,jj,kk)=(p(ii,jj,kk+1)-2*p(ii,jj,kk)+p(ii,jj,kk-1))/(dz.^2);
            
            advx_x(ii,jj,kk) = (advx(ii+1,jj,kk)-advx(ii-1,jj,kk))/2.0/dx;
            advy_y(ii,jj,kk) = (advy(ii,jj+1,kk)-advy(ii,jj-1,kk))/2.0/dy;
            advz_z(ii,jj,kk) = (advz(ii,jj,kk+1)-advz(ii,jj,kk-1))/2.0/dz;
        end
    end
end

advx_x(1,:,:) = (advx(2,:,:)-advx(1,:,:)) / dx;
advx_x(end,:,:) = (advx(end,:,:)-advx(end-1,:,:)) / dx;
advx_x(:,1,:) = (advx(:,2,:)-advx(:,1,:)) / dx;
advx_x(:,end,:) = (advx(:,end,:)-advx(:,end-1,:)) / dx;
advx_x(:,:,1) = (advx(:,:,2)-advx(:,:,1)) / dx;
advx_x(:,:,end) = (advx(:,:,end)-advx(:,:,end-1)) / dx;

advy_y(1,:,:) = (advy(2,:,:)-advy(1,:,:)) / dy;
advy_y(end,:,:) = (advy(end,:,:)-advy(end-1,:,:)) / dy;
advy_y(:,1,:) = (advy(:,2,:)-advy(:,1,:)) / dy;
advy_y(:,end,:) = (advy(:,end,:)-advy(:,end-1,:)) / dy;
advy_y(:,:,1) = (advy(:,:,2)-advy(:,:,1)) / dy;
advy_y(:,:,end) = (advy(:,:,end)-advy(:,:,end-1)) / dy;

advz_z(1,:,:) = (advz(2,:,:)-advz(1,:,:)) / dz;
advz_z(end,:,:) = (advz(end,:,:)-advz(end-1,:,:)) / dz;
advy_y(:,:,1) = (advy(:,:,2)-advy(:,:,1)) / dz;
advz_z(:,end,:) = (advz(:,end,:)-advz(:,end-1,:)) / dz;
advz_z(:,:,1) = (advz(:,:,2)-advz(:,:,1)) / dz;
advz_z(:,:,end) = (advz(:,:,end)-advz(:,:,end-1)) / dz;

b = -advx_x-advy_y-advz_z;
dxyz = dy*dy*dz*dz+dx*dx*dz*dz+dx*dx*dy*dy;
dxyz1 = dx*dx*dy*dy*dz*dz;

for ii=2:NX-1
    for jj=2:NY-1
        for kk=2:NZ-1            
            p1(ii,jj,kk) = -b(ii,jj,kk)*dxyz1/(2*dxyz);
        end
    end
end

% px(2:end-1,2:end-1,2:end-1) = (p(3:end,2:end-1,2:end-1) + 2*p(2:end-1,2:end-1,2:end-1) - p(1:end-2,2:end-1,2:end-1))/dx/dx; 
% px(1,2:end-1,2:end-1) = px(2,2:end-1,2:end-1); px(end,2:end-1,2:end-1) = px(end-1,2:end-1,2:end-1);
% 
% py(2:end-1,2:end-1,2:end-1) = (p(2:end-1,3:end,2:end-1) + 2*p(2:end-1,2:end-1,2:end-1) - p(2:end-1,1:end-2,2:end-1))/dy/dy; 
% py(2:end-1,1,2:end-1) = py(2:end-1,2,2:end-1); py(2:end-1,end,2:end-1) = py(2:end-1,end-1,2:end-1);
% 
% pz(2:end-1,2:end-1,2:end-1) = (p(2:end-1,2:end-1,3:end) + 2*p(2:end-1,2:end-1,2:end-1) - p(2:end-1,2:end-1,1:end-2))/dz/dz; 
% pz(2:end-1,2:end-1,1) = pz(2:end-1,2:end-1,2); pz(2:end-1,2:end-1,end) = pz(2:end-1,2:end-1,end-1);
% p1(2:end-1,2:end-1,2:end-1) = 


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
surf(advx(:,:,NZ/2))
% subplot(2,2,1); slice(x,y,z,advs_x,pi,pi,pi),colorbar
% subplot(2,2,2); slice(x,y,z,advs_y,pi,pi,pi),colorbar
% subplot(2,2,3); slice(x,y,z,advs_z,pi,pi,pi),colorbar
% subplot(2,2,4); slice(x,y,z,u,pi,pi,pi),colorbar

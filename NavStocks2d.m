clear all;
clc;
close all;

nu = 0.001;

NX = 128; NY = 128;
dx = 2*pi/NX; dy = 2*pi/NY;

dt = 1e-4; tTimes = 100000;

[x,y] = meshgrid(dx*[1:NX],dy*[1:NY]);
% 
% u = exp(-(x-pi).^2-(y-pi).^2);
% v = zeros(NX,NY);

a = 1; b = 1;
u = cos(a*x).*sin(b*y);
v = sin(a*x).*cos(b*y);

% u = random('unif',-1,1,NX,NY);
% v = random('unif',-1,1,NX,NY);

kx = ones(1,NY)'*[0:NX/2 -NX/2+1:-1];
ky = [0:NY/2 -NY/2+1:-1]'*ones(1,NX);

ikx = 1i*kx;
iky = 1i*ky;

ikx2 = ikx.*ikx;
iky2 = iky.*iky;

k2 = ikx.^2 + iky.^2;
k2p = k2;
k2p(1,1) = 1;

uf = fft2(u); vf = fft2(v);

for nt=1:tTimes  

    
    uxf = ikx.*uf; ux = real(ifft2(uxf));
    uyf = iky.*uf; uy = real(ifft2(uyf));
    uxxf = ikx2.*uf; 
    uyyf = iky2.*uf;
    
    vxf = ikx.*vf; vx = real(ifft2(vxf));
    vyf = iky.*vf; vy = real(ifft2(vyf));
    vxxf = ikx2.*vf; 
    vyyf = iky2.*vf;
    
    advxf = fft2(u.*ux+v.*uy);
    advyf = fft2(u.*vx+v.*vy);
    
    dissxf = -nu*k2.*uf;
    dissyf = -nu*k2.*vf;
    
    pf = -(ikx.*advxf+iky.*advyf)./k2p;
    p = real(ifft2(pf));
    
    uf_new = uf + dt*(-advxf-ikx.*pf+dissxf);
    vf_new = vf + dt*(-advyf-iky.*pf+dissyf);
    
    u = real(ifft2(uf_new));
    v = real(ifft2(vf_new));    
    
    if mod(nt,100) == 0        
        subplot(2,2,1); contourf(u,20),title(nt),colorbar,colormap('jet'),shading flat; drawnow
        subplot(2,2,2); contourf(v,20),title(nt*dt),colorbar,colormap('jet'),shading flat; drawnow
        subplot(2,2,3); contourf(p,20),title(nt*dt),colorbar,colormap('jet'),shading flat; drawnow
    end
    
    uf = uf_new; vf = vf_new;
end
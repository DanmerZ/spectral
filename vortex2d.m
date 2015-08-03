clear all;
clc;
close all;

nu = 0.001;

NX = 128; NY = 128;
dx = 2*pi/NX; dy = 2*pi/NY;

dt = 1e-3; tTimes = 30000;

[x,y] = meshgrid(dx*[1:NX],dy*[1:NY]);

w = exp(-((x-pi).^2+(y-pi+pi/4).^2)/(0.2))+exp(-((x-pi).^2+(y-pi-pi/4).^2)/(0.2))-0.5*exp(-((x-pi-pi/4).^2+(y-pi-pi/4).^2)/(0.4));
% w = cos(x) - cos(y);

kx = ones(1,NY)'*[0:NX/2 -NX/2+1:-1];
ky = [0:NY/2 -NY/2+1:-1]'*ones(1,NX);

ikx = 1i*kx;
iky = 1i*ky;

k2 = ikx.^2 + iky.^2;
k2p = k2;
k2p(1,1) = 1;

wf = fft2(w);
for nt=1:tTimes
    psif = -wf./k2p;
    vx = real(ifft2(iky.*psif));
    vy = -real(ifft2(ikx.*psif));
    
    wx = real(ifft2(ikx.*wf));
    wy = real(ifft2(iky.*wf));
    
    conv = vx.*wx + vy.*wy;
    convf = fft2(conv);
    
%     wf_new = wf + dt*(-convf+nu*k2.*wf);
    wf_new = (wf.*(1+.5*dt*nu*k2)-dt*convf) ./ (1-.5*dt*nu*k2);
    
    if mod(nt,25) == 0
        w = real(ifft2(wf_new));
        contourf(w,50),title(nt*dt),colorbar,colormap('jet'),shading flat; drawnow  
%         subplot(2,2,1); contourf(w,50),title(nt*dt),colorbar,colormap('jet'),shading flat; drawnow     
%         subplot(2,2,2); contourf(vx,50),title(nt*dt),colorbar,colormap('jet'),shading flat; drawnow
%         subplot(2,2,3); contourf(vy,50),title(nt*dt),colorbar,colormap('jet'),shading flat; drawnow
%         subplot(2,2,4); contourf(sqrt(vx.^2+vy.^2),50),title(nt*dt),colorbar,colormap('jet'),shading flat; drawnow
    end
    
    wf = wf_new;

end
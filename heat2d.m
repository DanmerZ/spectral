clear all;
clc;
close all;
alpha = 1;

NX = 128; NY = 128;
dx = 2*pi/NX; dy = 2*pi/NY;

[i,j] = meshgrid(1:NX,1:NY);
x = i*dx; y = j*dy;

t0 = 0.0; tSteps = 100000;
dt = 0.001*dx;

shfactor = 100;
nshow = tSteps / shfactor;

ikx = 1i*ones(1,NY)'*[0:NX/2 -NX/2+1:-1];
iky = 1i*[0:NY/2 -NY/2+1:-1]'*ones(1,NX);

ikx2 = ikx.*ikx;
iky2 = iky.*iky;

u = zeros(nshow,NX,NY);

uold = exp(-((x-pi).^2 + (y-pi).^2));
unew = uold;
u(1,:,:) = reshape(uold,1,NX,NY);



for tstep = 1:tSteps-1
    uold = unew;
    uf = fft2(uold); % fftn(squeeze(u(tstep,:,:)));
    uxx = real(ifft2(ikx2.*uf));
    uyy = real(ifft2(iky2.*uf));
    %u(tstep+1,:,:) = u(tstep,:,:) + ...
    %    dt*alpha*reshape(uxx.^2 + uyy.^2,1,NX,NY);     
    unew = uold + dt*alpha*(uxx + uyy);
    if (mod(tstep,shfactor) == 0 ) 
        u(tstep/shfactor,:,:) = unew;
        fprintf('%d\n',tstep)
    end    
end




for tstep=1:nshow-1
    surf(squeeze(u(tstep,:,:))),title(sprintf('t=%d',tstep)), axis([0 NX 0 NY 0 1]), drawnow;    
    M(tstep) = getframe;    
end





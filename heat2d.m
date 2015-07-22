clear all;
clc;
close all;
alpha = 1;

NX = 128; NY = 128;
dx = 2*pi/NX; dy = 2*pi/NY;

[i,j] = meshgrid(1:NX,1:NY);
x = i*dx; y = j*dy;

t0 = 0.0; tSteps = 150000;
dt = 0.001*dx;

shfactor = 300;
nshow = tSteps / shfactor;

ikx = 1i*ones(1,NY)'*[0:NX/2 -NX/2+1:-1];
iky = 1i*[0:NY/2 -NY/2+1:-1]'*ones(1,NX);

ikx2 = ikx.*ikx;
iky2 = iky.*iky;

u = zeros(nshow,NX,NY);
uu = zeros(nshow,NX,NY);

uold = exp(-((x-pi).^2 + (y-pi).^2));
unew = uold;
uuold = -exp(-((x-pi).^2 + (y-pi).^2)); % uold;
uunew = uuold;

u(1,:,:) = reshape(uold,1,NX,NY);
uu(1,:,:) = reshape(uuold,1,NX,NY);

for tstep = 1:tSteps-1
    uold = unew; uuold = uunew;
    uf = fftn(uuold); % fftn(squeeze(u(tstep,:,:)));
%     uxx = real(ifft2(ikx2.*uf));
%     uyy = real(ifft2(iky2.*uf));
%     [uxx,uyy] = sd2_2d(uf,ikx2,iky2);
 
    uunew = uuold + dt*rhs(uf,alpha,ikx2,iky2); % RK4(uuold, dt,alpha,ikx2,iky2); 
    unew = uold + dt*uunew;    
    
    if (mod(tstep,shfactor) == 0 ) 
        u(tstep/shfactor,:,:) = reshape(unew,1,NX,NY); 
        uu(tstep/shfactor,:,:) = reshape(uunew,1,NX,NY);
        fprintf('%d, %f\n',tstep,unew(1,1))
    end 
    
    % BC
     unew(:,1) = 0.0; unew(1,:) = 0.0;
     unew(:,end) = 0.0; unew(end,:) = 0.0;
end




for tstep=1:nshow-1
%     , axis([0 NX 0 NY 0 1]),
    surf(squeeze(u(tstep,:,:))),title(sprintf('t=%d',tstep)),axis([0 NX 0 NY 0 1]),  drawnow;    
    M(tstep) = getframe;    
end





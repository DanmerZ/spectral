clear all;
clc;
close all;
alpha = 1;

NX = 128; NY = 128;
dx = 2*pi/NX; dy = 2*pi/NY;

[i,j] = meshgrid(1:NX,1:NY);
x = i*dx; y = j*dy;

t0 = 0.0; tSteps = 1000;
dt = 0.005*dx;

shfactor = 100;
nshow = tSteps / shfactor;

ikx = 1i*ones(1,NY)'*[0:NX/2 -NX/2+1:-1];
iky = 1i*[0:NY/2 -NY/2+1:-1]'*ones(1,NX);

ikx2 = ikx.*ikx;
iky2 = iky.*iky;

u = zeros(nshow,NX,NY);

uu = zeros(tSteps,NX,NY);

uold = .1*exp(-100*((x-pi).^2 + (y-pi).^2));
unew = uold;
u(1,:,:) = reshape(uold,1,NX,NY);
uu(1,:,:) = reshape(uold,1,NX,NY);

for tstep = 1:tSteps-1
    uold = unew;
    uf = fft2(uold); % fftn(squeeze(u(tstep,:,:)));
%     uxx = real(ifft2(ikx2.*uf));
%     uyy = real(ifft2(iky2.*uf));
%     [uxx,uyy] = sd2_2d(uf,ikx2,iky2);
 
    unew = RK4(uold, dt,alpha,ikx2,iky2); %uold + dt*rhs(uf,alpha,ikx2,iky2);
    uu(tstep+1,:,:) = uu(tstep,:,:) + dt*reshape(unew,1,NX,NY);
    
    if (mod(tstep,shfactor) == 0 ) 
        u(tstep/shfactor,:,:) = unew;        
        fprintf('%d\n',tstep)
    end 
    
    % BC
%     unew(:,1) = 0.1; % unew(1,:) = 1.0;
%     unew(:,end) = -0.1; % unew(end,:) = -1.0;
end




for tstep=1:tSteps-1
%     , axis([0 NX 0 NY 0 1]),
    surf(squeeze(uu(tstep,:,:))),title(sprintf('t=%d',tstep)), axis([0 NX 0 NY -1 1]), drawnow;    
    M(tstep) = getframe;    
end





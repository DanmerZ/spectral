clear all;

nu = 1.0e-3; % viscosity
NX = 128;
NY = 128;
dt = 1e-1;
TF = 1000.0;
TSCREEN=25;

I = sqrt(-1);
dx = 2*pi/NX;
dy = 2*pi/NY;
t = 0.0;

% IC
w = random('unif',-1,1,NX,NY);

% matrix of wave numbers
kx = I*ones(1,NY)'*(mod((1:NX)-ceil(NX/2+1),NX)-floor(NX/2));
ky = I*(mod((1:NY)'-ceil(NY/2+1),NY)-floor(NY/2))*ones(1,NX);

dealias = kx<2/3*NX & ky<2/3*NY;

ksquare_viscous = kx.^2 + ky.^2;
ksquare_poisson = ksquare_viscous;
ksquare_poisson(1,1) = 1;

w_hat = fft2(w);

k = 0;
while t < TF
    k = k + 1;
    
    psi_hat = -w_hat./ksquare_poisson;
    
    u = real(ifft2( ky.*psi_hat ));
    v = real(ifft2( -kx.*psi_hat ));
    
    w_x = real(ifft2( kx.*w_hat ));
    w_y = real(ifft2( ky.*w_hat ));
    
    conv = u.*w_x + v.*w_y;
    conv_hat = fft2(conv);
    
    conv_hat = dealias.*conv_hat;
    
    w_hat_new = w_hat + dt*( nu*ksquare_viscous.*w_hat - conv_hat );
    
    t = t + dt;
    
    if k == TSCREEN
        w = real(ifft2(w_hat_new));
        contourf(w,50); colorbar; shading flat; colormap('jet');
        title(num2str(t));
        drawnow
        k = 0;
    end
    
    w_hat = w_hat_new;
    
end
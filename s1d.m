clear all;

%%matlab
N = 16;
x = linspace(2*pi/N,2*pi,N);
ik = 1i*fftshift([-N/2+1:N/2]); %1i*[0:N/2 -N/2+1:-1]; % i * wave number vector (matlab ordering)
ik2 = ik.*ik; % multiplication factor for second derivative
u = exp(cos(x));
u_hat = fft(u);
v_hat = ik2 .* u_hat;
v = real(ifft(v_hat)); 

plot(x,v- (sin(x).^2 - cos(x)) .* exp(cos(x)))


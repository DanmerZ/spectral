clear;

N = 24; 
x = cos(pi*(0:N)/N);
y = x';
dt = 6/N^2;
[xx,yy] = meshgrid(x,y);
plotgap = round((1/3)/dt);
dt = (1/3)/plotgap;

vv = exp(-40*((xx-.4).^2 + yy.^2));
vvold = vv;

% time-stepping leap-frog formula
[ay,ax] = meshgrid([.56 .06],[.1 .55]); clf
for n=0:3*plotgap
    t = n*dt;
    if rem(n+.5,plotgap) < 1
        i = n/plotgap+1;
        subplot('position',[ax(i) ay(i) .36 .36])
        [xxx,yyy] = meshgrid(-1:1/16:1,-1:1/16:1);
        vvv = interp2(xx,yy,vv,xxx,yyy,'cubic');
        mesh(xxx,yyy,vvv), axis([-1 1 -1 1 -0.15 1])
        colormap([0 0 0]), title(num2str(t)), drawnow
    end
    uxx = zeros(N+1,N+1); uyy = zeros(N+1,N+1);
    ii = 2:N;
    for i = 2:N
        v = vv(i,:); V = [v fliplr(v(ii))];
        U = real(fft(V));
        W1 = real(ifft(1i*[0:N-1 0 1-N:-1].*U));
        W2 = real(ifft(-[0:N 1-N:-1].^2.*U));
        uxx(i,ii) = W2(ii)./(1-x(ii).^2) - x(ii).*W1(ii)./(1-x(ii).^2).^(3/2);
    end
    for j = 2:N
        v = vv(:,j); V = [v; flipud(v(ii))];
        U = real(fft(V));
        W1 = real(ifft(1i*[0:N-1 0 1-N:-1]'.*U));
        W2 = real(ifft(-[0:N 1-N:-1]'.^2.*U));
        uyy(ii,j) = W2(ii)./(1-y(ii).^2) - y(ii).*W1(ii)./(1-y(ii).^2).^(3/2);
    end
    vvnew = 2*vv - vvold + dt^2*(uxx+uyy);
    vvold = vv; vv = vvnew;
end

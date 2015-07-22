function [ unew ] = RK4( uold, dt,alpha,ikx2,iky2  )
%RK4 4-th order Runge-Kutta in time

uf = fftn(uold);
k1 = rhs(uf,alpha,ikx2,iky2);

uf = fftn(uold + .5*dt*k1);
k2 = rhs(uf,alpha,ikx2,iky2);

uf = fftn(uold + .5*dt*k2);
k3 = rhs(uf,alpha,ikx2,iky2);

uf = fftn(uold + dt*k3);
k4 = rhs(uf,alpha,ikx2,iky2);

unew = uold + (1.0/6.0)*dt*(k1+2*k2+2*k3+k4);


end


function [ u,t ] = euler1( f,t0,tf,u0,n )
%[u,t]=euler1(f,t0,tf,u0,n)
%
% This function implements Euler's method for solving the IVP
%
% du/dt=f(t,u), u(t0)=u0
%
% on the interval [t0,tf]. n steps of Euler's method are taken;
% the step size is dt=(tf-t0)/n.
% Compute the grid and allocate space for the solution
k = length(u0);
t=linspace(t0,tf,n+1);
u=zeros(k,n+1);

u(:,1) = u0;
dt = (tf-t0)/n;

for ii=1:n
   u(:,ii+1)=u(:,ii)+dt*f(t(ii),u(:,ii)); 
end

end


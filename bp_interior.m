function out=bp_interior(D,y,maxiter)
% BPINTERIOR - Basis Pursuit demostration
% 
% Reference :
%    S.Chen, Atomic Decomposition by Basis Pursuit, SIAM, 2001
% 
% Copyright (C) Liangfu Chen (2012-2013)

% clear all;clf;
if nargin<1,D = [1/2 sqrt(3)/2; 1 0; -1/sqrt(2) -1/sqrt(2)]';end
if nargin<2,y = [1 1/2]';else,y=y';end
if nargin<3,maxiter=2;end
assert(size(D,1)==size(y,1));
assert(1==size(y,2));

% clear all;clf;
% D = [1/2 sqrt(3)/2; 1 0; -1/sqrt(2) -1/sqrt(2)]';
% y = [1 1/2]';
% maxiter=2;

A=D;%[D -D];
b=y;
c=ones([size(A,2) 1]);
z=ones([size(A,2) 1])*.1;
x=ones([size(A,2) 1])*.1;
mu=ones([size(A,2) 1])*.1;
e=ones([size(A,2) 1]);
gamma=1e-4;
delta=1e-4;
I=eye(size(A,2));

for ii=1:maxiter
%% (a) compute residuals and diagonal matrix D
t=c+(gamma^2)*x-z-A'*y;
r=b-A*x-(delta^2)*y;
v=mu.*e-diag(z)*x;
D=inv(inv(diag(x))*diag(z)+(gamma^2)*I);

%% (b) solve least square approximation
dy=(A*D*A'+(delta^2)*eye(length(b)))\(r+A*D*(t-inv(diag(x))*v));
dx=D*(A'*dy+inv(diag(x))*v-t);
dz=inv(diag(x))*(v-diag(z)*dx);

%% (c) update the variables
stepx=1e20; stepz=1e20;
blocking=find(dx<0);
if length(blocking)>0,stepx=min(x(blocking)./(-dx(blocking)));end
blocking=find(dz<0);
if length(blocking)>0,stepz=min(z(blocking)./(-dz(blocking)));end
rho_p=min(.99*stepx,1);
rho_d=min(.99*stepz,1);
x=x+rho_p*dx;
y=y+rho_d*dy;
z=z+rho_d*dz;
mu=(1-min([rho_p,rho_d,.99]))*mu;

%% termination criteria
if ii==2,break;end
end

out.val=x;
out.loc=1:length(x);

% A,x,A*x

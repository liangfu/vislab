function out=OrthogonalMatchingPursuit(D,y,maxiter,epsilon)
% ORTHOGONALMATCHINGPURSUIT - Orthogonal Matching Pursuit demo
% Copyright (C) Liangfu Chen (2012-2013)

% clear all;
if nargin<1,D = [1/2 sqrt(3)/2; 1 0; -1/sqrt(2) -1/sqrt(2)]';end
if nargin<2,y = [1 1/2]';end
if nargin<3,maxiter=10;end
if nargin<4,epsilon=0.02;end

y=y(:);
assert(size(D,1)==size(y,1));
phi=[];
%% estimation
err=y;
for iter=1:maxiter
[maxval,maxloc]=max(err'*D); % inner-product
out.loc(iter)=maxloc;
phi=[phi D(:,maxloc)];
x=inv(phi'*phi+1e-8)*(phi'*y); % least square
out.val=x;
err=y-phi*x;
if(sqrt(sum(err.^2))/size(y,1)<epsilon),break;end
end

%% evaluation by reconstructing the signal
% approx=D(:,out.loc(:))*out.val(:);
% [approx,y]

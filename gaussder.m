function gaussder(X,Y,sigma)
% GAUSSDER computes gaussian derivatives
% 
% Description:
%    gaussder(X,Y,SIGMA)
% 
% See also:
%    gauss, gabor2

if nargin<3,sigma=2;end
if nargin<1,sz=16;[X,Y]=meshgrid(-(sz-1)/2.:(sz-1)/2.);end

G=exp(-(X.^2+Y.^2)./(2*sigma^2));
dx=imfilter(G,[-1,0,1]);
dy=imfilter(G,[-1;0;1]);
subplot(2,5,1),imshow(dx,[]);
subplot(2,5,2),imshow(dy,[]);

dxdx=imfilter(dx,[-1,0,1]);
subplot(2,5,3),imshow(dxdx,[]);

dxdx=imfilter(dx,[-1,0,0;0,0,0;0,0,1]);
subplot(2,5,4),imshow(dxdx,[]);

dxdx=imfilter(dx,[0,0,-1;0,0,0;1,0,0]);
subplot(2,5,5),imshow(dxdx,[]);

subplot(2,5,6),imshow(G,[]);

dxdx=imfilter(dx,[-1,0,1]);
dxdxdx=imfilter(dxdx,[-1,0,1]);
subplot(2,5,7),imshow(dxdxdx,[]);

dxdx=imfilter(dx,[-1,0,0;0,0,0;0,0,1]);
dxdxdx=imfilter(dxdx,[-1,0,0;0,0,0;0,0,1]);
subplot(2,5,8),imshow(dxdxdx,[]);

dxdx=imfilter(dx,[0,0,-1;0,0,0;1,0,0]);
dxdxdx=imfilter(dxdx,[0,0,-1;0,0,0;1,0,0]);
subplot(2,5,9),imshow(dxdxdx,[]);

dydy=imfilter(dy,[-1;0;1]);
dydydy=imfilter(dydy,[-1;0;1]);
subplot(2,5,10),imshow(dydydy,[]);

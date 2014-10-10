function g = gauss(X,Y,sigma)
% GAUSS generates 2D gaussion kernel 
% 
% Example:
%    [X,Y]=meshgrid(-10:10,-10:10);
%    sigma=15;
%    g=guass(X,Y,sigma);
%    pcolor(X,Y,g);

if nargin==0,[X,Y]=meshgrid(-10:10,-10:10);sigma=15;end

g=exp(-0.5*(X.^2+Y.^2)./sigma);

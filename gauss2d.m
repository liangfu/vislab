function G=gauss2d(X,Y,theta,sigma,center)
% GAUSS2D generate 2d rotated gaussian kernel
% Example:
%    [X,Y]=meshgrid(-3:.1:3,-3:.1:3);theta=pi/6;sigma=[2,1];
%    g=gauss2d(X,Y,theta,sigma);
sigma_x = sigma(1);
sigma_y = sigma(2);
if nargin<5,x0 = 0; y0 = 0;else,x0=center(1);y0=center(2);end
a = cos(theta)^2/2/sigma_x^2 + sin(theta)^2/2/sigma_y^2;
b = -sin(2*theta)/4/sigma_x^2 + sin(2*theta)/4/sigma_y^2 ;
c = sin(theta)^2/2/sigma_x^2 + cos(theta)^2/2/sigma_y^2;

G = exp( - (a*(X-x0).^2 + 2*b*(X-x0).*(Y-y0) + c*(Y-y0).^2)) ;

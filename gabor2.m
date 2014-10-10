function gb=gabor2(x_range,y_range,sigma,theta,psi,gamma)
% GABOR2 create 2D gabor kernel
% gb=gabor2(x_range,y_range,sigma,theta,psi,gamma)

if nargin<6,if nargin<5,
if nargin<4,if nargin<3,
if nargin<2,if nargin<1,
x_range=-5:.1:5;end;
y_range=x_range;end;
sigma=1;end;
theta=pi/3;end;
psi=0;end;
gamma=1;end

sigma_x = sigma;
sigma_y = sigma/gamma;
 
% Bounding box
[x,y] = meshgrid(x_range,y_range);
 
% Rotation 
x_theta=x*cos(theta)+y*sin(theta);
y_theta=-x*sin(theta)+y*cos(theta);
 
gb= exp(-.5*(x_theta.^2/sigma_x^2+y_theta.^2/sigma_y^2)).*...
    cos(2*pi*x_theta+psi);

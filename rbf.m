function kappa = rbf(x_i,x_j,epsilon,theta1, theta2, theta3, beta)
% RBF - generate radial basis function
if nargin<4, theta1=1.0;theta2=1.0;theta3=0.0;beta=1.0; end
invbeta=beta^-1;
tmp = theta1.*exp(-0.5*theta2*((x_i-x_j).^2))+theta3+invbeta*epsilon;
kappa = sum(sum(tmp));

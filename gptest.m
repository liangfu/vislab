function gptest

xs = (-5:0.2:5)'; ns = size(xs,1); 
phi(1)=1;phi(2)=1;phi(3)=0;phi(4)=1e9;
fs = m(xs) + chol(rbf(xs,xs,phi))'*randn(ns,1);
plot(xs,fs,'.'),colormap('jet');

function hyp=gptrain(x,y)

function mu=m(x) 
  mu=0;%.25*x.^2;

function kappa=rbf(p, q, phi)
  if nargin<3, theta1=1.0;theta2=1.0;theta3=0.0;beta=1.0; end
  tmp = repmat(p',size(q))-repmat(q, size(p'));
  invbeta = phi(4)^-1;
  kappa = phi(1).*exp(-0.5*phi(2)*tmp.^2)+phi(3)+invbeta*eye(size(p,1));

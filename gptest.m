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

  function K = covfnGE(x,y,par)
% COVFNGE Gamma Exponential Covariance Function
%K = covfnGE(x,y,par)
K = par(1).value.*exp(-par(2).value*(sqdist(x,y)).^(0.5*par(3).value));
if length(par)==4
    K = K + par(4).value*eye(size(K,1));
end

function l=logdet(A)
% LOGDET Log determinant of a positive definite matrix computed in a
%numerically more stable manner
[u s v]=svd(A); 
l=sum(log(diag(s)+1.0e-20));

function D = sqdist(x,y)
% SQDIST Square distance between vectors in x and y
% D = sqdist(x,y)
[d m]=size(x); n=size(y,2);
xx=repmat(sum(x.^2,1),n,1); 
yy=repmat(sum(y.^2,1),m,1); 
xy=x'*y; 
D = xx'+yy-2*xy;
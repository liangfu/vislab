function [meantrain,meantest,logpygx]=GPclass(xtrain,ctrain,xtest)
% GPCLASS Gaussian Process Binary Classification
% 
% [meantrain,meantest,logpygx]=GPclass(xtrain,ctrain,xtest)
% 
% See also demoGPclass

% define the Kernel matrix:
par(1).value=1; % prefactor
par(2).value=2; % inverse root lengthscale
par(3).value=2; % gamma exponential
par(4).value=0.001; % jitter
opts.maxits = 100; opts.tol = 1e-8;

N=size(xtrain,2);
train=1:N; test=(1:size(xtest,2))+N;
% Covariance function for noisy outputs
K = feval('covfnGE',[xtrain xtest],[xtrain xtest],par); 

ctrain=ctrain(:); % make a column vector
y=zeros(N,1);
for wloop=1:opts.maxits % Newton update for Laplace approximation
    sigvec=sigma(y);
    D = diag(sigvec.*(1-sigvec));
    tmp = D*y+ ctrain-sigvec;
    Mat=eye(N)+K(train,train)*D;
    yold = y;
    y=K(train,train)*(Mat\tmp);
    if mean(abs(y-yold))<opts.tol; break; end
end
meantrain=sigvec;

% compute predictions:
cms = ctrain-sigvec; invD = diag(1./diag(D));
for n=1:size(xtest,2)
    k = K(train,test(n));
    kss = K(test(n),test(n));
    mn = k'*cms;  % mean of the field projection
    vr = kss - k'*((K(train,train)+invD)\k); % variance
    meantest(1,n)=avsigmaGauss(mn,vr);
end

% training data log likelihood
logpygx = ctrain'*y - sum(log(1+exp(y))) - 0.5*y'*cms-0.5*logdet(Mat);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                           UTILITY                            %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s=sigma(x)
% SIGMA 1./(1+exp(-x))
% s=sigma(x) = 1./(1+exp(-x))
s=1./(1+exp(-x));

function out=avsigmaGauss(mn,v)
% AVSIGMAGAUSS Average of a logistic sigmoid under a Gaussian
% out=avsigmaGauss(mn,v)
% simple approximation for the average of 1/(1+exp(-x)) over a Gaussian
% with mean mn and variance v
erflambda=sqrt(pi)/4;
out=0.5+0.5*erf(erflambda*mn/sqrt(1+2*erflambda^2*v));

function K = covfnGE(x,y,par)
% COVFNGE Gamma Exponential Covariance Function
%K = covfnGE(x,y,par)
K = par(1).value.*exp(-par(2).value*(sqdist(x,y)).^(0.5*par(3).value));
if length(par)==4
    K = K + par(4).value*eye(size(K,1));
end

function D = sqdist(x,y)
% SQDIST Square distance between vectors in x and y
% D = sqdist(x,y)
[d m]=size(x); n=size(y,2);
xx=repmat(sum(x.^2,1),n,1); 
yy=repmat(sum(y.^2,1),m,1); 
xy=x'*y; 
D = xx'+yy-2*xy;

function l=logdet(A)
% LOGDET Log determinant of a positive definite matrix computed in a
%numerically more stable manner
[u s v]=svd(A); 
l=sum(log(diag(s)+1.0e-20));


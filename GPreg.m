function [meantest,vartest,logpygx]=GPreg(xtrain,ytrain,xtest,par,covfn)
%GPREG Gaussian Process Regression
% [meantest,vartest,logpygx]=GPreg(xtrain,ytrain,xtest,par,covfn)
%
% Inputs:
% xtrain : matrix of training data. Each column contains a datapoint
% ytrain : vector of 1 dimensional outputs corresponding to xtrain
% xtest  : the test inputs
% par : covariance parameters
% covfn : covariance function
% 
% See demoGPreg.m, gppred.m
N=size(xtrain,2);
% Covariance function for noisy outputs
K = feval('covfnGE',xtrain,xtrain,par); 
for t=1:size(xtest,2)
	xs=xtest(:,t);
	Kxs=feval('covfnGE',xtrain,xs,par(1:3));
	Kss=feval('covfnGE',xs,xs,par(1:3));
	tmp = K\Kxs;
	meantest(:,t) = ytrain(:)'*tmp;
	vartest(:,t) = Kss - tmp'*Kxs;
end

% training data log likelihood
logpygx = -0.5*ytrain(:)'*K/ytrain(:)'-0.5*logdet(K)-0.5*N*log(2*pi);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                           UTILITY                            %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

function demoGPreg
%DEMOGPREG demo of Gaussian Process regression
figure

par(1).value=1; % prior variance of the clean output
par(2).value=2; % inverse scale
par(3).value=2; % gamma for the Gamma-Exponential Covariance
par(4).value=0.0001; % jitter
covfn='covfnGE'; % Gamma-Exponential Covariance Function

% try to understand the GP function prior by drawing sample functions
x=-10:0.05:10; P=size(x,2);
K = feval(covfn,x,x,par); % Covariance function
y=mvrandn(zeros(P,1),K,5);
plot(x,y);

% make some training data:
xtrain=[(-0.75 + rand(1,20)) (0.75+ rand(1,20))]; N=size(xtrain,2);
ytrain = sin(4*xtrain)+ 0.1*randn(1,N);

xtest=[-4:0.1:4];
[meantest,vartest,logpygx]=GPreg(xtrain,ytrain,xtest,par,covfn);

figure
subplot(221),plot(xtest,meantest,'r-');  title('mean prediction');
subplot(222),plot(xtest,meantest-sqrt(vartest),'g-'); title('standard error bars')
subplot(223),plot(xtest,meantest+sqrt(vartest),'g-');
subplot(224),plot(xtrain,ytrain,'.');

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
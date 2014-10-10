function gpc

clear all; clf;
X=[];Y=[];
if 1,
N = [40,20]; % number of samples
D = 2; % dimension of feature vector
sz = [N(1),1];
tmp=randn(sz)+5;
angle=rand(sz)*2*pi-pi;
X=[tmp.*cos(angle), tmp.*sin(angle);randn([N(2),D])]+...
  repmat((randn([1,2])*1.2+10),[sum(sum(N)),1]);
Y=[zeros([N(1),1]);ones([N(2),1])];
idx=randperm(sum(N));
X=X(idx,:);Y=Y(idx,:);
end

N=size(X,1);
maxiter=2;
sigma=inline('1./exp(-0.5*x)');
K = inline('exp(-0.5*(repmat(p'',size(q))-repmat(q,size(p''))).^2)');

f=zeros([N,1]);
for iter=1:maxiter
y=sigma(f);
logpy=-log(y);
W=-(y==1).*(1-(y==1));
L=chol(eye(N)+sqrt(W)'*K(y',y)*sqrt(W));
b=W.*f+((y+1)/2.-(y==1));
a=b-(sqrt(W)'*L')\(L\(sqrt(W)'*K*b));
f_old=f;
f=K*a;
break;
if norm(f-f_old)<1e-2,break;end
end

% logqy=-.5*a'*f+logpy-sum(sum(log(L)));


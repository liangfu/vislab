function PRIOR=GPClearn(X,Y)
% GPCLEARN learning gaussian process classification

if nargin<1
X=[13.14,10.85;11.74,6.83;17.42,8.73;18.32,13.57;8.31,7.60;
   8.69,8.83;8.78,6.03;16.19,15.13;12.37,11.90;11.69,6.34;
   13.78,11.93;13.33,11.92;13.44,7.89;12.30,12.32;9.68,12.40];
Y=[1;0;0;0;0;0;0;0;1;0;1;1;0;1;0];
end

cases=unique(Y);
w1=find(Y==cases(1));w2=find(Y==cases(2));
Y(w1)=-1;Y(w2)=1;cases=[-1,1];
maxiter=length(Y);
N=size(X,1);
K=exp(-0.5*(sum(X.^2,2)*ones(1,N)+ones(N,1)*sum(X.^2,2)'-2*(X*X')));

f=zeros(N,1);sumlogL=0;
for iter=1:maxiter
PI=1./(1+exp(-f(:,end))); 
W=diag(PI.*(1-PI));
L=chol(eye(N)+sqrt(W)*K*sqrt(W)');
b=W*f(:,end)+((Y+1)/2-PI);
a=b-sqrt(W)*L'\(L\(sqrt(W)*K*b));
f(:,iter)=K*a;
disp(sprintf('%d:%f',iter,norm(f(:,iter)-f(:,max(1,iter-1)))))
if norm(f(:,iter)-f(:,max(1,iter-1)))<1e-5 && iter>1,break;end
end
f=f(:,end);
disp(sprintf('converged with %d iterations',iter));

PRIOR.X=X;
PRIOR.Y=Y;
PRIOR.alpha=Y-PI;
PRIOR.f = f;
PRIOR.logpygx=-.5*a'*f+PI-.5*log(det(eye(N)+sqrt(W)*K*sqrt(W)'));

kern=inline(['exp(-0.5*(sum(X0.^2,2)*ones(1,size(X1,1))+' ...
             'ones(size(X0,1),1)*sum(X1.^2,2)''-2*(X0*X1'')))'],'X0','X1');
L = real(kern(X,PRIOR.X)*PRIOR.alpha);
[mean(Y(find(L<0))) mean(Y(find(L>0)))]

%% display
w1=find(L<0);
w2=find(L>0);
subplot(121),plot(X(w1,1),X(w1,2),'r*',X(w2,1),X(w2,2),'b*');
subplot(122),hist(L,100);%plot(L(w1,1),L(w1,1),'r.',L(w2,1),L(w2,1),'b.');

function PRIOR=SGKFDlearn(X,Y)
% SGKFDLEARN learning Sparse Greedy Kernel Fisher Discriminant
% 
% implemented accorrding to
% Mika, Sebastian, A. J. Smola, and Bernhard Schölkopf. "An
% improved training algorithm for kernel fisher discriminants."
% Proceedings AISTATS. Vol. 2001. Morgan Kaufmann, 2001.
% 
% Copyright (C) Liangfu Chen (2012-2013)
    
warning off;
opts1= optimset('display','off');

if nargin<1,
X=[13.14,10.85;11.74,6.83;17.42,8.73;18.32,13.57;8.31,7.60;
   8.69,8.83;8.78,6.03;16.19,15.13;12.37,11.90;11.69,6.34;
   13.78,11.93;13.33,11.92;13.44,7.89;12.30,12.32;9.68,12.40];
Y=[1;0;0;0;0;0;0;0;1;0;1;1;0;1;0]; 
end
maxiter=length(Y);

cases=unique(Y);
assert(2==length(cases)); % number of classes
ell=size(X,1);

kern=inline(['exp(-0.5*(sum(X0.^2,2)*ones(1,size(X1,1))+' ...
             'ones(size(X0,1),1)*sum(X1.^2,2)''-2*(X0*X1'')))'],'X0','X1');
K=exp(-0.5*(sum(X.^2,2)*ones(1,ell)+ones(ell,1)*sum(X.^2,2)'-2*(X*X')));

w1=find(Y==cases(1));
w2=find(Y==cases(2)); 
Y(w1)=1;Y(w2)=-1;cases(1)=1;cases(2)=-1;
N1=length(w1);
N2=length(w2);
one_mat=ones([ell,1]);
one_pos=max(Y,0);
one_neg=max(-Y,0);

m=0;

%% compute between and within class scatters
I=[];alpha=[];err=[];
ellset=1:ell;
for iter=1:maxiter
S=ellset(~ismember(1:ell,I));
objmax=inf;
for ii=S % numerate elements
I_=[I,ii];%I_=union(I,ii);
K_=K(:,I_);
ell_=size(K_,2); 
H_ = [ell          one_mat'*K_
      K_'*one_mat  K_'*K_+1e-5*eye(ell_)]; 
      % K_'*one_mat  K_'*K_+eye(ell_)*1e-2];
%% update inverse of H
if (size(H_,2)<3)
invH_=inv(H_);
else
B=H_(1:end-1,end);C=H_(end,end);gamma=inv(C-B'*invH*B);invHB=invH*B;
invH_=[invH+invHB*gamma*invHB', -gamma*invHB;
       -(gamma*invHB)',          gamma];
end
A_p = [N1; K_'*one_pos];
A_n = [N2; K_'*one_neg];
c=[N1-N2; K_'*Y];
%% dual variables
% lambda=quadprog([A_p'*invH_*A_p A_p'*invH_*A_n 
%                  A_n'*invH_*A_p A_n'*invH_*A_n],...
%                 -[-N1+c'*invH_*A_p;N2+c'*invH_*A_n],...
%                 [],[],[],[],[],[],[],opts1);
lambda=inv([A_p'*invH_*A_p A_p'*invH_*A_n 
            A_n'*invH_*A_p A_n'*invH_*A_n])*...
       [-N1+c'*invH_*A_p;N2+c'*invH_*A_n];
lambda_p=lambda(1);lambda_n=lambda(2);
%% primal variables
% a=quadprog(H_,-c,[],[],[A_p';A_n'],[N1;-N2],[],[],[],opts1);
a=invH_*(c-(lambda_p*A_p+lambda_n*A_n));
%% dual objective
dualobj = -.5*a'*H_*a+[-N1;N2]'*[lambda_p;lambda_n]+.5*ell;
if dualobj<objmax,ii_opt=ii;objmax=dualobj;end
end

I=[I,ii_opt];%union(I,ii_opt);
K_=K(:,I);
H = [ell          one_mat'*K_
     K_'*one_mat  K_'*K_+1e-5*eye(size(K_,2))]; 
     % K_'*one_mat K_'*K_+eye(size(K_,2))*1e-2];
%% update inverse of H
if (size(H,2)<3)
invH=inv(H);
else
B=H(1:end-1,end);C=H(end,end);gamma=inv(C-B'*invH*B);invHB=invH*B;
invH=[invH+invHB*gamma*invHB', -gamma*invHB;
       -(gamma*invHB)',          gamma];
end

% a=quadprog(H,-c,[],[],[A_p';A_n'],[N1;-N2],[],[],[]);
a=invH*(c-(lambda_p*A_p+lambda_n*A_n));
b=a(1);  % threshold
alpha=a(2:end); % coefficient
%% termination criteria
% err(iter)=sqrt(sum(sum((K_*alpha+one_mat*b-Y).^2)))/length(Y);
objmaxarr(iter)=objmax;
err(iter)=abs(objmaxarr(end)-objmaxarr(max(1,end-1)));
disp(sprintf('%d:%f',iter,err(end)));
% L = real(kern(X,X(I,:))*alpha)+b;
% disp(sprintf('%d:%f at [%.3f,%.3f] with b=%f',...
%              iter,err(end),mean(Y(find(L<0))),mean(Y(find(L>0))),b));
% sumerr(iter)=sum([mean(Y(find(L<0))),mean(Y(find(L>0)))]);
if mean(err(max(1,end-int32(maxiter/20)):end))<.02 && iter>1,break;end
end

PRIOR.X=X(I,:);
PRIOR.alpha=alpha;
PRIOR.b=b;

%% project to lower dimension
L = real(kern(X,PRIOR.X)*alpha)+b;
[mean(Y(find(L<0))) mean(Y(find(L>0)))]

%% display
w1=find(L<0);
w2=find(L>0);
subplot(121),plot(X(w1,1),X(w1,2),'r*',X(w2,1),X(w2,2),'b*');
subplot(122),hist(L,100);%plot(L(w1,1),L(w1,1),'r.',L(w2,1),L(w2,1),'b.');



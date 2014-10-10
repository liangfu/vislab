function PRIOR=KFDlearn(X,Y)
% KFDLEARN learning Kernel Fisher Discriminant
    
if nargin<1,
X=[13.14,10.85;11.74,6.83;17.42,8.73;18.32,13.57;8.31,7.60;
   8.69,8.83;8.78,6.03;16.19,15.13;12.37,11.90;11.69,6.34;
   13.78,11.93;13.33,11.92;13.44,7.89;12.30,12.32;9.68,12.40];
Y=[1;0;0;0;0;0;0;0;1;0;1;1;0;1;0];
end

C=unique(Y);
nclasses = length(C); % number of classes
sumN=size(X,1);

K=exp(-0.5*(sum(X.^2,2)*ones(1,sumN)+ones(sumN,1)*sum(X.^2,2)'-2*(X*X')));

w1=find(Y==C(1));
w2=find(Y==C(2)); 
N1=length(w1);
N2=length(w2);
one_mat=ones([sumN,1]);
mu1 = ((1/N1)*K*(Y==C(1)))';
mu2 = ((1/N2)*K*(Y==C(2)))';

%% compute between and within class scatters
M = (mu2-mu1)*(mu2-mu1)';
N = K'*K-(N1*mu1*mu1'+N2*mu2*mu2');
% regularization
N = N+eye(size(N))*1e-8*trace(N);

%% compute alpha
[eigvec,eigval]=eig(inv(N)*(mu2-mu1)'*(mu2-mu1));
[~,idx]=sort(real(diag(eigval)),'descend');
eigvec=eigvec(:,idx);eigval=eigval(:,idx);
alpha=eigvec(:,1:2);

%% project to lower dimension
L = real(K*alpha);

PRIOR.alpha=alpha;
PRIOR.X=X;

%% display
subplot(121),plot(X(w1,1),X(w1,2),'r*',X(w2,1),X(w2,2),'b*');
subplot(122),plot(L(w1,1),L(w1,2),'r.',L(w2,1),L(w2,2),'b.');


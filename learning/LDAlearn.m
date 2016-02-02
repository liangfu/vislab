function [W,Prior]=lda(X,Y,D_low)

%% generate data sets
if nargin<2
N = [120,80,40]; % number of samples
K = length(N); % number of classes
D = 50; % dimension of feature vector
sigma = rand([K,D]).*1.0;
X=[];Y=[];
for i=1:K
  X = [X;randn(N(i),D)+repmat(sigma(i,:),[N(i),1])];  
  Y = [Y;ones([N(i),1])*(i-1)];
end
idx=randperm(size(X,1));
X=X(idx,:);Y=Y(idx,:);
end

C=unique(Y);
K=length(C);
D=size(X,2);
sumN=size(X,1);
if nargin<3,D_low=min(2,D/2);end

%% compute mean,covariance etc.
clear dat;
for i=1:K
  dat(i).w=find(Y==C(i));
  dat(i).mu=mean(X(dat(i).w,:));
  dat(i).S= ...
           (X(dat(i).w,:)-repmat(dat(i).mu,[length(dat(i).w) 1]))'* ...
           (X(dat(i).w,:)-repmat(dat(i).mu,[length(dat(i).w) 1]));
  dat(i).N=length(dat(i).w);
end

%% mean of all classes
mu=mean(X,1);

%% within-class scatter 
Sw=zeros(size(dat(1).S));
for i=1:K,Sw = Sw+dat(i).S; end
Sw=Sw+eye(size(Sw))*1e-8;

%% between-class scatter
Sb=zeros([D D]);
for i=1:K,Sb=Sb+dat(i).N*(dat(i).mu-mu)'*(dat(i).mu-mu);end

%% eigen decomposition
[U S V]=svd(inv(Sw)*Sb);%size(U)
W=U(:,1:D_low)';
Prior=zeros(K,D_low);
for i=1:K, Prior(i,:)=dat(i).mu*W'; end

%% compute probability
L = X*W';
P = zeros(size(L,1),2);
for i=1:K
  tmp = exp(-0.5*sum((L-repmat(Prior(i,:),[sumN,1])).^2,2));
  P(:,i) = tmp(:,1);
end
P = P./repmat(sum(P,2),[1,size(P,2)]);

%% display
if 1
clf;
separate_ratio = [mean(P(Y==C(1),:));mean(P(Y==C(2),:))]
subplot(121),
% imagesc(real(separate_ratio));%,colormap('gray');
plot(X(Y==C(1),1),X(Y==C(1),2),'r*',...
     X(Y==C(2),1),X(Y==C(2),2),'b*'),title('input space');
subplot(122),
% plot(X(Y==C(1),1),X(Y==C(1),2),'r*',...
%      X(Y==C(2),1),X(Y==C(2),2),'b*');
plot(L(Y==C(1),1),L(Y==C(1),2),'r*',...
     L(Y==C(2),1),L(Y==C(2),2),'b*'),title('feature space');
end

% if nargout==0,W=[];end

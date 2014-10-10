function lda3

%% generate data sets
N = [10,10,10]; % number of samples
K = length(N); % number of classes
D = 2; % dimension of feature vector
sigma = rand([K,D]).*15.0;
X=[];Y=[];
for i=1:K
  X = [X;randn(N(i),D)+repmat(sigma(i,:),[N(i),1])];  
  Y = [Y;ones([N(i),1])*(i-1)];
end
C=unique(Y);

%% compute mean,covariance etc.
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

%% between-class scatter
Sb=zeros([D D]);
for i=1:K,Sb=Sb+dat(i).N*(dat(i).mu-mu)'*(dat(i).mu-mu);end

%% eigen decomposition
[U S V]=svd(inv(Sw)*Sb);
W=U(:,1:2)';
Prior=zeros(K,2);
for i=1:K, Prior(i,:)=dat(i).mu*W'; end

%% compute probability
L = X*W';
for i=1:K
    P(:,i) = exp(-0.5*sum((L-repmat(Prior(i,:),[sum(N),1])).^2,2));
end
P = P./repmat(sum(P,2),[1,3]);

%% display
fres = [];
for i=1:K
  fres=[fres;mean(P(sum(N(1:i-1))+1:sum(N(1:i)),:))];
end
fres
subplot(121),
imagesc(fres),colormap('gray');
subplot(122),
xmat=[];
for i=1:K,
  xmat=[xmat, ...
        X(sum(N(1:i-1))+1:sum(N(1:i-1))+N(i),1)+ ...
        X(1+sum(N(1:i-1)):sum(N(1:i-1))+N(i),2)*j];
end
plot(xmat,'.');


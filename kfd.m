function opt=kfd(X,Y)
% KFD - Kernel Fisher Discriminant Analysis

%% generate data sets
if nargin~=2,
X=[];Y=[];
N = [10,5]; % number of samples
D = 2; % dimension of feature vector
if 0
sigma = rand([length(N),D]).*5.0+8.0;
X=[];Y=[];
for i=1:length(N)
  X = [X;randn(N(i),D)+repmat(sigma(i,:),[N(i),1])];  
  Y = [Y;ones([N(i),1])*(i-1)];
end
else
sz = [N(1),1];
tmp=randn(sz)+5;
angle=rand(sz)*2*pi-pi;
X=[tmp.*cos(angle), tmp.*sin(angle);randn([N(2),D])]+...
  repmat((randn([1,2])*1.2+10),[sum(sum(N)),1]);
Y=[zeros([N(1),1]);ones([N(2),1])];
end
idx=randperm(sum(N));
X=X(idx,:);Y=Y(idx,:);
end

%% size of training data
C=unique(Y);
nclasses = length(C); % number of classes
sumN=size(X,1);

%% construct using Gaussian kernel 
% K = zeros([sumN,sumN]);
% for r = 1:sumN
% for c = 1:r
%   K(r,c) = exp(-0.5*sum(sum(((X(r,:) - X(c,:)).^2))));
%   K(c,r) = K(r,c);
% end
% end
K=exp(-0.5*(sum(X.^2,2)*ones(1,sumN)+ones(sumN,1)*sum(X.^2,2)'-2*(X*X')));

%% compute mean,covariance etc.
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
N = N+eye(size(N))*1e-3*trace(N);

%% compute alpha
[eigvec,eigval]=eig(inv(N)*(mu2-mu1)'*(mu2-mu1));
[~,idx]=sort(real(diag(eigval)),'descend');
eigvec=eigvec(:,idx);eigval=eigval(:,idx);
alpha=eigvec(:,1:2);

%% project to lower dimension
L = real(K*alpha);

if 1
%% testing
% generate grid data for testing
minX=min(X);
maxX=max(X);
bound=(maxX-minX)*0.1;
range_row=linspace(minX(1)-bound(1),maxX(1)+bound(1),30);
range_col=linspace(minX(2)-bound(2),maxX(2)+bound(2),30);
[X_map_test,Y_map_test]=meshgrid(range_row,range_col);
sumN_test=size(X_map_test,1)*size(X_map_test,2);
X_test = [reshape(X_map_test,[sumN_test,1]), ...
          reshape(Y_map_test,[sumN_test,1]),zeros([sumN_test,size(X,2)-2])];
% start testing 
% K_test = zeros([sumN_test,sumN]);
% for r = 1:sumN_test
% for c = 1:sumN
%   K_test(r,c) = exp(-0.5*sum(sum(((X_test(r,:) - X(c,:)).^2))));
% end
% end
K_test=exp(-0.5*(sum(X_test.^2,2)*ones(1,sumN)+...
                 ones(sumN_test,1)*sum(X.^2,2)'-2*(X_test*X')));
proj_test=K_test*alpha;
proj_test(proj_test(:,1)>max(max(L(:,1))),1)=max(max(L(:,1)));
proj_test(proj_test(:,1)<min(min(L(:,1))),1)=min(min(L(:,1)));

%% display 
clf;
subplot(121)
hold on,
colormap('gray');
axis([minX(1)-bound(1),maxX(1)+bound(1),minX(2)-bound(2),maxX(2)+bound(2)]),
pcolor(range_row,range_col,reshape(proj_test(:,1),size(X_map_test))),
plot(X(w1,1),X(w1,2),'r*',X(w2,1),X(w2,2),'b*');
hold off;
subplot(122)
end
if 0
  hist(L(w1,1));hold on;hist(L(w2,1));
  h = findobj(gca,'Type','patch');% display(h);
  set(h(1),'FaceColor','r','EdgeColor','k');
  set(h(2),'FaceColor','b','EdgeColor','k');
else
  plot(L(w1,1),L(w1,2),'r.',L(w2,1),L(w2,2),'b.');
end

opt.alpha=alpha;
opt.L=L;
opt.X=X;

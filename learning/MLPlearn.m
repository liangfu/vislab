function net = MLPlearn(net,X,Y)
% MLPLEARN - MLP For Binary Classification

% generate training data sets
if nargin==1
N = [145,135]; % number of samples
K = length(N);   % number of classes
D = 3;          % dimension of feature vector
sigma = rand([K,D]).*1.0+eye([K,D]).*1.0;
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
maxiter=net.maxiter; % maximum iteration
hsize=round((D+K)*.5) % number of neurons in hidden layer
net.layer{1}.w=rand([hsize,D])*.1;
if length(net.layer)>2
  for l=2:length(net.layer)-1,net.layer{l}.w=rand([hsize,hsize])*.1;end
end
net.layer{end}.w=rand([K,hsize])*.1;
net.layer{1}.input = X';
net.layer{end}.target = [Y';1-Y'];

errors=[];
for iter=1:maxiter
% forward pass
for l=1:length(net.layer)
  net.layer{l}.output = net.layer{l}.forward(net.layer{l});
  if l<length(net.layer),net.layer{l+1}.input = net.layer{l}.output;end
end
% backward pass - step 1
net.layer{end}.dE=(net.layer{end}.target-net.layer{end}.output);
for l=length(net.layer):-1:1
  net.layer{l}.delta = net.layer{l}.backward_step1(net.layer{l});
  if l>1,net.layer{l-1}.dE = (net.layer{l}.w'*net.layer{l}.delta);end
end
% backward pass - step 2
for l=1:length(net.layer)
  net.layer{l}.w = net.layer{l}.backward_step2(net.layer{l});
end
%net.layer{l}
errors(iter)=abs(sum(sum(net.layer{end}.delta)));
end
subplot(121),plot(errors),title('error');

% generate test data set
if nargin==0
X=[];Y=[];
for i=1:K
  X = [X;randn(N(i),D)+repmat(sigma(i,:),[N(i),1])];  
  Y = [Y;ones([N(i),1])*(i-1)];
end
idx=randperm(size(X,1));
X=X(idx,:);Y=Y(idx,:);
end
% foward pass
net.layer{1}.input = X';
for l=1:length(net.layer)
  net.layer{l}.output = net.layer{l}.forward(net.layer{l});
  if l<length(net.layer),net.layer{l+1}.input = net.layer{l}.output;end
end
% labeling
L1=find((net.layer{end}.output(1,:))'<(net.layer{end}.output(2,:))'); sizeL1 = size(L1)
L2=find((net.layer{end}.output(1,:))'>(net.layer{end}.output(2,:))'); sizeL2 = size(L2)
subplot(122),plot(X(L1,1),X(L1,2),'r.',X(L2,1),X(L2,2),'b.');

end


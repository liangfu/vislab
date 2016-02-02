function net = MLPlearn(X,Y)
% RNNLEARN - RNN For Binary Classification

% generate training data sets
if nargin==0
N = [120,80]; % number of samples
K = length(N);   % number of classes
D = 5;          % dimension of feature vector
sigma = rand([K,D]).*20.0;
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
hsize=round((D+K)*.5); % number of neurons in hidden layer
alpha=.001;  % learning rate
maxiter=200; % maximum iteration

% hyperbolic
net.layer{1}.w=rand([hsize,D])*.1;
net.layer{1}.bias=0;
% logistic
net.layer{2}.w=rand([2,hsize])*.1;
net.layer{2}.bias=0;

errors=[];
for iter=1:maxiter
net.layer{1}.forward = @(x) tanh(net.layer{1}.w*x+net.layer{1}.bias);
net.layer{2}.forward = @(x) 1./(1+exp(-(net.layer{2}.w*x+net.layer{2}.bias)));
% forward pass
h = net.layer{1}.forward(X');
b = net.layer{2}.forward(h); 
% backward pass
net.layer{2}.delta = b.*(1-b).*([Y';1-Y']-b);
net.layer{1}.delta = (1-h.*h).*(net.layer{2}.w'*net.layer{2}.delta);
net.layer{2}.w = net.layer{2}.w + alpha * (net.layer{2}.delta * h');
net.layer{1}.w = net.layer{1}.w + alpha * (net.layer{1}.delta * X);
errors(iter)=abs(sum(sum(net.layer{2}.delta)));
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
h = net.layer{1}.forward(X');
b = net.layer{2}.forward(h);
% labeling
L1=find((b(1,:))'<(b(2,:))'); sizeL1 = size(L1)
L2=find((b(1,:))'>(b(2,:))'); sizeL2 = size(L2)
subplot(122),plot(X(L1,1),X(L1,2),'r.',X(L2,1),X(L2,2),'b.');

end


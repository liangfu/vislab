function net = RNNlearn(X,Y)
% RNNLEARN - RNN for sequence classification

% generate training data sets
if nargin==0
sz = [1,1380];I = zeros(sz);
I(110:220)=rand([1,111])*.1+1;
X=real(ifft(I));X=X(250:end-251).*1000;
X=([X;sin(3:2+size(X,2))])';
Y=zeros(size(X));Y(:,1)=1;
% subplot(121),plot(X,'-')
% subplot(122),plot(abs(fft(X)),'-')
end

C=unique(Y);
K=length(C);
D=30;%size(X,2);
hsize=round((D+K)*.5); % number of neurons in hidden layer
alpha=.001;  % learning rate
maxiter=200; % maximum iteration

pts = unique(randi(size(X,1)-D-1,[1,200]));
samples = [];
for ii=1:length(pts)
  
end

% hyperbolic
net.layer{1}.Wih=rand([hsize,D])*.1; % for current input 
net.layer{1}.Whh=rand([hsize,D])*.1; % for previous activation
net.layer{1}.Whk=rand([hsize,hsize])*.1; % for output
% logistic
net.layer{2}.w=rand([2,hsize])*.1;

errors=[];
for iter=1:maxiter
% forward pass
b0 = zeros([D,1]);
b1 = tanh(net.layer{1}.Wih*x+net.layer{1}.hh*b0);
b2 = 1./(1+exp(-(net.layer{2}.w*x+net.layer{2}.bias)));
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


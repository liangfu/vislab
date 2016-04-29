function MLPdemo()

net.alpha=0.002;
net.maxiter=200;

% hyperbolic - input layer
net.layer{1}.forward = @tanh_forward;
net.layer{1}.backward = @tanh_backward;
net.layer{1}.alpha=net.alpha;

% logistic - output layer
net.layer{2}.forward = @sigmoid_forward;
net.layer{2}.backward = @sigmoid_backward;
net.layer{2}.alpha=net.alpha;

% softmax - output layer for classification
% net.layer{2}.forward = @softmax_forward;
% net.layer{2}.backward = @softmax_backward;
% net.layer{2}.alpha=net.alpha;

% generate training data sets
N = [45,35,20]; % number of samples
K = length(N);  % number of classes
D = 4;          % dimension of feature vector
sigma = rand([K,D]).*1.0+eye([K,D]).*5.0;
X=[];Y=[];
for i=1:K
  X = [X;randn(N(i),D)+repmat(sigma(i,:),[N(i),1])];  
  Y = [Y;ones([N(i),1])*(i-1)];
end
idx=randperm(size(X,1)); % shuffle data list
X=X(idx,:);Y=Y(idx,:);
cidx=Y+1+(1:size(Y,1))'*K-K;
Z=zeros([K,size(cidx,1)]);Z(cidx)=1;
X=X';Y=Z;
tmp=X;X={};X{1}=tmp;clear tmp;

% training parameters
C=1:size(Y,1);
K=length(C);
D=size(X{1},1);
maxiter=net.maxiter; % maximum iteration
hsize=6; % number of neurons in hidden layer
net.layer{1}.W=rand([hsize,D])*.1; % initial weights
net.layer{2}.W=rand([K,hsize])*.1;

disp(sprintf('hidden layer size: %d',hsize));
disp(sprintf('number of classes: %d',K));

loss=[];
for iter=1:maxiter
[net.layer{1},X{2}]=net.layer{1}.forward(net.layer{1},X{1});

% gradchecks
% tanh_gradcheck(X{2},Y);
% sigmoid_gradcheck(X{2},Y);
% softmax_gradcheck(X{2},Y);

[net.layer{2},X{3}]=net.layer{2}.forward(net.layer{2},X{2}); 
dE_dY{3} = Y-X{3};
[net.layer{2},dE_dY{2}]=net.layer{2}.backward(net.layer{2},X{2},dE_dY{3});
[net.layer{1},dE_dY{1}]=net.layer{1}.backward(net.layer{1},X{1},dE_dY{2});
loss(iter)=sum(sum((dE_dY{3}.*dE_dY{3})));
end

if length(loss)>0,plot(loss);res=[(X{3})' Y'];res(1:20,:),end

end

function [layer,Y]=tanh_forward(layer,X)
layer.WX=layer.W*X;
Y=tanh(layer.WX);
end

function [layer,dE_dX]=tanh_backward(layer,X,dE_dY)
Y=tanh(layer.WX);
dE_dY_afder=(1-Y.*Y).*dE_dY;
dE_dX=dE_dY_afder'*layer.W;
layer.dE_dW=dE_dY_afder*X'; % for gradient checking
layer.W=layer.W+layer.alpha*layer.dE_dW;
end

function y=sigmoid(x) 
y=1./(1+exp(-x));end

function [layer,Y]=sigmoid_forward(layer,X)
layer.WX=layer.W*X;
Y=sigmoid(layer.WX);
end

function [layer,dE_dX]=sigmoid_backward(layer,X,dE_dY)
Y=sigmoid(layer.WX);
dE_dY_afder=Y.*(1-Y).*dE_dY;
dE_dX=(dE_dY_afder'*layer.W)';
layer.dE_dW=dE_dY_afder*X';
layer.W=layer.W+layer.alpha*layer.dE_dW;
end

function y=softmax(x)
expX=exp(x);
sumexpX=repmat(sum(expX),3,1);
y=expX./sumexpX;
end

function [layer,Y]=softmax_forward(layer,X)
layer.WX=layer.W*X;
Y=softmax(layer.WX);
end

function [layer,dE_dX]=softmax_backward(layer,X,dE_dY)
Y=softmax(layer.WX);
dE_dY_afder=Y.*(1-Y).*dE_dY;
dE_dX=(dE_dY_afder'*layer.W)';
layer.dE_dW=dE_dY_afder*X';
layer.W=layer.W+layer.alpha*layer.dE_dW;
end

function tanh_gradcheck(X,target)
X=X(:,1);target=target(:,1);
layer.W=rand([size(target,1),size(X,1)])*.1;
layer.forward = @tanh_forward;
layer.backward = @tanh_backward;
layer.alpha=0.01;
err=gradcheck(layer,X,target);
disp(sprintf('   tanh gradcheck error: %g',err));
end

function sigmoid_gradcheck(X,target)
X=X(:,1);target=target(:,1);
layer.W=rand([size(target,1),size(X,1)])*.1;
layer.forward = @sigmoid_forward;
layer.backward = @sigmoid_backward;
layer.alpha=0.01;
err=gradcheck(layer,X,target);
disp(sprintf('sigmoid gradcheck error: %g',err));
end

function softmax_gradcheck(X,target)
X=X(:,1);target=target(:,1);
layer.W=rand([size(target,1),size(X,1)])*.1;
layer.forward = @softmax_forward;
layer.backward = @softmax_backward;
layer.alpha=0.01;
err=gradcheck(layer,X,target);
disp(sprintf('softmax gradcheck error: %g',err));
end

function error=gradcheck(layer,X,target)
delta=1e-4;
W=layer.W;
[layer,Y]=layer.forward(layer,X);
[layer,dE_dX]=layer.backward(layer,X,target-Y);
loss=zeros(size(W));
for ridx=1:size(W,1),for cidx=1:size(W,2)
layer_p.W=W;W(ridx,cidx)=W(ridx,cidx)+delta;
[layer_p,Y]=layer.forward(layer_p,X);loss_p=sum(sum((target-Y).^2));
layer_m.W=W;W(ridx,cidx)=W(ridx,cidx)-delta;
[layer_m,Y]=layer.forward(layer_m,X);loss_m=sum(sum((target-Y).^2));
loss(ridx,cidx)=(loss_p-loss_m)/(delta*2);
end;end
error=mean(reshape(abs(loss-layer.dE_dW),1,[]));
end



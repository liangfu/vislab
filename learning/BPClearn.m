function model=BPClearn(X,Y)
% BPC data classification using backpropagation
% model=bpc(X,Y)
% 
% input:
%  - X     input data
%  - Y     label 

% generate data sets
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

classreg=0;
train0=[X zeros(sum(N),K)];
for ii=1:K,idx=find(Y==(ii-1));train0(idx,D+ii)=1;end
[Kr,MN]=size(train0);
% train0(1:10,:) % first 10 rows of training set
M=D,N=K
K=min(64,Kr)
t = 1;  % initialize epoch counter
ptr=1;  % initialize pointer for re-sampling the training file
not_converged = 1;  % not yet converged
scalein=1;
bstrate=0;      % initialize classification rate on tuning set to 0
nck=1;
nstall=0;   % initilize # of no improvement count. when nstall >
            % maxstall, quit.
maxstall=10;
alphaval=0.1; % learning rate (between 0 and 1, default = 0.1)
mom=0.8; % momentum constant (between 0 and 1, default 0.8)
L1=2;
L=L1+1; 
disp(['The output layer has ' int2str(N) ' neurons.']);
n(1) = M;
n(L) = N; % # of output neurons = output dimension
% w is a cell array with i-th element being a weight matrix
% of the i-th layer neurons.
if L > 2, 
   for i=2:L-1,
      n(i)=40;%input(['Number of neurons in hidden layer # ' int2str(i-1) ' = ']);
      w{i}=0.001*randn(n(i),n(i-1)+1); % first column is the bias weight
      dw{i}=zeros(size(w{i}));         % initialize dw
   end
end
w{L}=0.005*randn(n(L),n(L-1)+1); % first column is the bias weight
dw{L}=zeros(size(w{L}));
% ==============================================================
% choose types of activation function
% default: hidden layers, tanh (type = 2), output layer, sigmoid (type = 1)
% default parameter T = 1 is used.
% ==============================================================
atype=2*ones(L,1);  atype(L)=1; % default
chostype=0;
% ==============================================================
% termination criteria
% A. Terminate when the max. # of epochs to run is reached.
% ==============================================================
nepoch=100;
tune=train0;
chos=1;
% '1 - Use the entire training set to estimate training error;          '
% '2 - Use a separate fixed tuning data file to estimate training error;'
% '3 - Partition training set dynamically into training and tuning sets;'
% '    (This is for pattern classification problem)                     '];
test0=train0;
labeled=1;

% BP iterations begins
while not_converged==1,
	% start a new epoch
   % Randomly select K training samples from the training set.
   [train,ptr,train0]=rsample(train0,K,Kr,ptr); % train is K by M+N
   z{1}=(train(:,1:M))';   % input sample matrix  M by K
   d=train(:,M+1:MN)';     % corresponding target value  N by K
   
   % Feed-forward phase, compute sum of square errors 
   for l=2:L,              % the l-th layer
      u{l}=w{l}*[ones(1,K);z{l-1}]; % u{l} is n(l) by K
      z{l}=actfun(u{l},atype(l));
   end
   error=d-z{L};           % error is N by K
   E(t)=sum(sum(error.*error));
   
   % Error back-propagation phase, compute delta error 
   delta{L}=actfunp(u{L},atype(L)).*error;  % N (=n(L)) by K
   if L>2,
      for l=L-1:-1:2,
         delta{l}=(w{l+1}(:,2:n(l)+1))'*delta{l+1}.*actfunp(u{l},atype(l));
      end
   end
   
   % update the weight matrix using gradient, momentum and random perturbation
   for l=2:L,
      dw{l}=alphaval*delta{l}*[ones(1,K);z{l-1}]'+mom*dw{l}+randn(size(w{l}))*0.005;
      w{l}=w{l}+dw{l};
   end

   % display the training error
   %bpdisplay;

   % Test convergence to see if the convergence condition is satisfied, 
   % cvgtest;
% ==============================================================
% termination criteria
% A. Terminate when the max. # of epochs to run is reached.
% ==============================================================
if t>=nepoch
   not_converged=0;
   disp('Terminate training since Max # of epochs has been reached');
end
% ==============================================================
% B. Check the tuning set testing result periodically. If the tuning set testing
%    results is reducing, save the weights. When the tuning set testing results
%    start increasing, stop training, and use the previously saved weights.
% ==============================================================
if rem(t,nck)==0, % time to check tuning set testing result
   [Cmat,crate]=bptest(w,tune,atype);%disp(Cmat);
   disp(['Epoch ' int2str(t) ' classification rate : ' num2str(crate)]);
   if crate > bstrate, % if the classification of tuning set is improving
       bstrate=crate;
       wbest=w; % memorize the best weights
       nstall=0; % reset the no-improvement count
   else
       nstall=nstall+1; % if no improvement, increment the no-improvement count
   end % if crate
   if nstall>maxstall, % continuous maxstall no improvement, terminate
      not_converged=0;
      disp(['Terminate because no improvement for ' int2str(nstall) ' consecutive checks']);
   end % if nstall
   % for pattern classificaiton problem, if 
   % tuning set classification rate = 100%, terminate
   if classreg==0 & crate == 100,
       not_converged=0;
       disp(['Terminate because classification rate of tuning set is 100%'])
   end
end


	 t = t + 1;    % increment epoch count
end % while loop

disp('Final training results:')
[Cmat,Crate]=bptest(wbest,tune,atype);disp('Cmat');
disp(Cmat);disp('Crate');disp(Crate);

disp('Apply trained MLP network to the testing data. The results are: ');
[Cmat,crate,cout]=bptest(wbest,test0,atype,labeled,N);
if labeled==1, 
  disp('Confusion matrix Cmat = '); disp(Cmat);
  disp(['classification = ' num2str(crate) '%']);
end


end % end of main function



function [Cmat,crate,cout]=bptest(w,test,atype,labeled,N)
% Usage: [Cmat,crate,cout]=bptest(w,test,atype,labeled,N)
% This routine is used for pattern classification problems only.
% testing MLP network for classification rate and confusion matrix
% test (K by M+N, if labeled=1, K by M if labeled = 0) 
%      the feature vectors (K by M) and if labeled=1, target vectors (K by N)
% w{l}, 1<= l <= L consists of the L-1 weight matrices, w{1}=[].
%       each w{l} is n(l) x (n(l-1)+1)
% assume the output uses 1 of N encoding. Thus, class #1 is [1 0 ... 0]
%       class #2 is [0 1 0 .. 0], etc.
% atype(1:L): specify activation function types for each layer
%       if ignored, type(l) = 1
% labeled=1 (default) if test data set has labels, = 0 otherwise
% N: label dimension, will not be needed if labeled =1 
% Cmat: N by N the confusion matrix
% crate = sum(diag(Cmat))/K: classification rate
% cout: classifiers output wrt each testing data, K by N
% used by cvgtest.m
%
% copyright (c) 1996-2001 by Yu Hen Hu
% Last modified: 3/20/2001, 10/15/2003

if nargin<4, labeled=1; end % default the test set has labels
if nargin<3, atype=ones(length(w),1); end

[K,MN]=size(test); % find # of samples and input+output dimension
[m1,M1]=size(w{2});  % find input dimension
M=M1-1; 
if labeled ==1, N=MN-M; end
L=length(w);

z{1}=(test(:,1:M))';     % input sample matrix  M by K
if labeled==1, 
   target=test(:,M+1:MN)';  % corresponding target value  N by K
end
   
% Feed-forward phase,  
for l=2:L,               % the l-th layer
  z{l}=actfun(w{l}*[ones(1,K);z{l-1}],atype(l)); % z{l} is n(l) by K
end

% To compute confusion matrix, we need at least two outputs. If the MLP has
% only one output to encode two classes, we change it to a 1 out of 2 encoding
if N==1,
   cout=[[z{L}>0.5];[z{L}<0.5]]+zeros(2,K);  % convert into one of N encoding
   if labeled==1, % if testing set has labels
      target=[target;ones(1,K)-target]; 
   end
   N=2;  % change size of confusion matrix to 2 by 2
else
   cout = [(z{L} - ones(N,1)*max(z{L})) == 0]+zeros(N,K); % N x K
end

if labeled==1, 
   % Then we calculate the NxN confusion matrix.
   Cmat=([(target-ones(N,1)*max(target))==0]+zeros(N,K))...
      *cout'; 
   crate=sum(diag(Cmat))*100/K;
elseif labeled==0, % if the testing set has no labels, only provide cout
   Cmat=eye(N); crate=1;
end

end % end of function bptest

function [train,ptr,x]=rsample(x,K,Kr,ptr)
% Usage: [train,ptr]=rsample(x,K,Kr,ptr)
% random sample K out of Kr rows from matrix x. 
% ptr: current row count pointer. K rows will be taken from
%      row # ptr+1 to ptr + K
%      if ptr + K > Kr, the x matrix will be randomized by row
%      and then the first K rows will be taken.
%      Each time ptr will be updated
% call: randomize.m
% (C) copyright 2001 by Yu Hen Hu
% created: 3/8/2001
if ptr+K > Kr,
   x=randomize(x);
   train=x(1:K,:);
   ptr=K;
else
   train=x(ptr+1:ptr+K,:);
   ptr=ptr+K;
end
end

function y=actfun(x,type,par)
% Usage: y=actfun(x,type,par)
% Compute activation functions
% x: net function, a K x N matrix
% y: activation function, a K x N matrix
% type: 1 - sigmoid, 2 - tanh, 3- linear, 4 - radial
% par: parameter list
%    sigmoid: T,  y=1/(1+exp(-x/T)), yp=y*(1-y)/T
%    tanh: T,     y=(exp(x/T)-exp(-x/T))/(exp(x/T)+exp(-x/T))
%                 yp=(1-y*y)/T
%    linear:m,b   y=ax+b, yp=a
%    radial:m,sig y=exp(-(x-m)^2/sig^2), yp=-2(x-m)*y/sig^2
% (C) copyright 2001 by YU HEN HU
% created: 2/5/2001
if nargin<=2,% if omitted from input vari list, set default
   if type==3,
      par=[1 0]; 
   elseif type==4,
      par=[0 1];
   else
      par=1;
   end
end
if nargin==1, % no type info either
   type=1;
end
switch type
case 1 % sigmoid
   T=par(1);
   y = 1./(1+exp(-x/T));
case 2 % tanh'
   T=par(1);
   tmp=exp(x/T);
   y=(tmp-1./tmp)./(tmp+1./tmp);
case 3 % linear'
   a=par(1); b=par(2);
   y=a*x+b;
case 4 % radial'
   m=par(1); sig=par(2);
   s=sig^2;
   tmp=(x-m).*(x-m);
   y=exp(-tmp/s);
end
end % end of function definition

function yp=actfunp(x,type,par)
% Usage: yp=actfunp(x,'type',par)
% Compute activation functions and their derivatives
% x: net function, a K x N matrix
% y: activation function, a K x N matrix
% type: 1 - sigmoid, 2 - tanh, 3- linear, 4 - radial
% par: parameter list
%    sigmoid: T,  y=1/(1+exp(-x/T)), yp=y*(1-y)/T
%    tanh: T,     y=(exp(x/T)-exp(-x/T))/(exp(x/T)+exp(-x/T))
%                 yp=(1-y*y)/T
%    linear:m,b   y=ax+b, yp=a
%    radial:m,sig y=exp(-(x-m)^2/sig^2), yp=-2(x-m)*y/sig^2
% (C) copyright 2001 by YU HEN HU
% created: 2/5/2001
if nargin<=2,% if omitted from input vari list, set default
   if type==3,
      par=[1 0]; 
   elseif type==4,
      par=[0 1];
   else
      par=1;
   end
end
if nargin==1, % no type info either
   type=1;
end
switch type
case 1 % sigmoid
   T=par(1);
   y = 1./(1+exp(-x/T));
   yp = y.*(ones(size(y))-y)/T;
case 2 % tanh'
   T=par(1);
   tmp=exp(x/T);
   y=(tmp-1./tmp)./(tmp+1./tmp);
   yp=(ones(size(y))-y.*y)/T;
case 3 % linear'
   a=par(1); b=par(2);
   y=a*x+b;
   yp=a*ones(size(y));
case 4 % radial'
   m=par(1); sig=par(2);
   s=sig^2;
   tmp=(x-m).*(x-m);
   y=exp(-tmp/s);
   yp=(-2/sig^2)*(x-m).*y;
end
end % end of function definition

function [B,I]=randomize(A,rowcol)
% Usage: [B,I]=randomize(A,rowcol)
% randomize row orders or column orders of A matrix
% rowcol: if =0 or omitted, row order (default)
%         if = 1, column order
% copyright (C) 1996-2001 by Yu Hen Hu
% Last modified: 8/30/2001
rand('state',sum(100*clock))
if nargin == 1,
   rowcol=0;
end
if rowcol==0, 
   [m,n]=size(A);
   p=rand(m,1);
   [p1,I]=sort(p);
   B=A(I,:);
elseif rowcol==1,
   Ap=A';
   [m,n]=size(Ap);
   p=rand(m,1);
   [p1,I]=sort(p);
   B=Ap(I,:)';
end
end % end of function definition


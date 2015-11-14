function model=bpc(X,Y)
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
      dw{l}=alphaval*delta{l}*[ones(1,K);z{l-1}]'+...
            mom*dw{l}+randn(size(w{l}))*0.005;
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
   disp(['Epoch # ' int2str(t)]);
   [Cmat,crate]=bptest(w,tune,atype);disp(Cmat);
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




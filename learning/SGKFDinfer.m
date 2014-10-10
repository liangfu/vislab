function SGKFDinfer(PRIOR,X_test)
% SGKFDINFER inference with Sparse Greedy Kernel Fisher Discriminant

if nargin==0,
PRIOR=SGKFDlearn;
%% import learned prior info
alpha=PRIOR.alpha;
X=PRIOR.X;
sumN=size(X,1);

%% generate testing sample points
X_=[13.14,10.85;11.74,6.83;17.42,8.73;18.32,13.57;8.31,7.60;
   8.69,8.83;8.78,6.03;16.19,15.13;12.37,11.90;11.69,6.34;
   13.78,11.93;13.33,11.92;13.44,7.89;12.30,12.32;9.68,12.40];
minX=min(X_);
maxX=max(X_);
bound=(maxX-minX)*0.1;
range_row=linspace(minX(1)-bound(1),maxX(1)+bound(1),30);
range_col=linspace(minX(2)-bound(2),maxX(2)+bound(2),30);
[X_map_test,Y_map_test]=meshgrid(range_row,range_col);
sumN_test=size(X_map_test,1)*size(X_map_test,2);
X_test = [reshape(X_map_test,[sumN_test,1]), ...
          reshape(Y_map_test,[sumN_test,1]),zeros([sumN_test,size(X,2)-2])];
elseif nargin==2,
%% import learned prior info
alpha=PRIOR.alpha;
X=PRIOR.X;
sumN=size(X,1);
sumN_test=size(X_test,1);
else assert(false);
end

%% inference with testing points and prior knowledge
K_test=exp(-0.5*(sum(X_test.^2,2)*ones(1,sumN)+...
                 ones(sumN_test,1)*sum(X.^2,2)'-2*(X_test*X')));
L_test=K_test*alpha+PRIOR.b;

%% display
% for tt=-.5:.1:1.2
tt=0;
w1=find(L_test(:,1)>tt);
w2=find(L_test(:,1)<tt);
subplot(121),plot(X_test(w1,1),X_test(w1,2),'r*',...
                  X_test(w2,1),X_test(w2,2),'b*');
subplot(122),hist(L_test,100);
% subplot(122),plot(L_test(w1,1),L_test(w1,1),'r.',...
%                   L_test(w2,1),L_test(w2,1),'b.');
% drawnow,pause(5)
% end


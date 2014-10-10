function KFDinfer(PRIOR,X_test)
% KFDINFER inference with Kernel Fisher Discriminant

if nargin==0,
PRIOR=KFDlearn;
%% import learned prior info
alpha=PRIOR.alpha;
X=PRIOR.X;
sumN=size(X,1);

%% generate testing sample points
minX=min(X);
maxX=max(X);
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
kern=inline(['exp(-0.5*(sum(X0.^2,2)*ones(1,size(X1,1))+' ...
             'ones(size(X0,1),1)*sum(X1.^2,2)''-2*(X0*X1'')))'],'X0','X1');
% K_test=exp(-0.5*(sum(X_test.^2,2)*ones(1,sumN)+...
%                  ones(sumN_test,1)*sum(X.^2,2)'-2*(X_test*X')));
K_test=kern(PRIOR.X,X_test);
L_test=K_test'*alpha;

%% display
% for tt=-.5:.1:1.2
tt=0;
w1=find(L_test(:,1)>tt);
w2=find(L_test(:,1)<tt);
subplot(121),plot(X_test(w1,1),X_test(w1,2),'r*',...
                  X_test(w2,1),X_test(w2,2),'b*');
subplot(122),plot(L_test(w1,1),L_test(w1,2),'r.',...
                  L_test(w2,1),L_test(w2,2),'b.');
% drawnow,pause(5)
% end


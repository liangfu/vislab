function GPCinfer
% GPCINFER inference with gaussian process prior

PRIOR=GPClearn_old;

X=PRIOR.X;
alpha=PRIOR.alpha;
% invD=PRIOR.invD;

N=size(X,1);
minX=min(X);
maxX=max(X);
bound=(maxX-minX)*0.1;
range_row=linspace(minX(1)-bound(1),maxX(1)+bound(1),30);
range_col=linspace(minX(2)-bound(2),maxX(2)+bound(2),30);
[X_map_test,Y_map_test]=meshgrid(range_row,range_col);
sumN_test=size(X_map_test,1)*size(X_map_test,2);
X_test = [reshape(X_map_test,[sumN_test,1]), ...
          reshape(Y_map_test,[sumN_test,1]),zeros([sumN_test,size(X,2)-2])];
N_test=size(X_test,1);

%% inference with testing points and prior knowledge
% K=exp(-0.5*(sum(X.^2,2)*ones(1,N)+...
%                  ones(N,1)*sum(X.^2,2)'-2*(X*X')));
K_test=exp(-0.5*(sum(X_test.^2,2)*ones(1,N)+...
                 ones(N_test,1)*sum(X.^2,2)'-2*(X_test*X')));
% Kss_test=exp(-0.5*(sum(X_test.^2,2)*ones(1,N_test)+...
%                    ones(N_test,1)*sum(X_test.^2,2)'-2*(X_test*X_test')));
mn=K_test*alpha;

% vr=Kss_test'-K_test*((K+invD)\K_test');
% erflambda=sqrt(pi)/4;meantest=[];
% for ii=1:N_test
% meantest(ii,:)=...
% 2.*erf((erflambda*mn(ii))/sqrt(1+(2*erflambda^2)*vr(:,ii)));
% end
% [diag(meantest),alpha(6:end),Y(6:end)]

%% display
% for tt=-1.5:.01:1.2
L_test=mn;
tt=0;
w1=find(L_test(:,1)>tt);
w2=find(L_test(:,1)<tt);
subplot(121),plot(X_test(w1,1),X_test(w1,2),'r*',...
                  X_test(w2,1),X_test(w2,2),'b*');
subplot(122),plot(L_test(w1,1),L_test(w1,1),'r.',...
                  L_test(w2,1),L_test(w2,1),'b.');
% drawnow;
% end



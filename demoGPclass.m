function demoGPclass
%DEMOGPCLASS:
% demo of Gaussian Process Binary Classification using Laplace with fixed
% hyperparameters. Note that this is not optimally encoded for speed nor accuracy

% make some training data
M=20;
xtrain=randn(2,M); % class 1
xtrain=[xtrain (randn(2,M)+2.5*ones(2,M))]; % class 0
xtrain=[xtrain (randn(2,M)-2.5*ones(2,M))]; % class 0
ctrain(M+1:3*M)=1; ctrain(1:M)=0; % training class labels
N=3*M; % total number of training points

% define predictions for on a very coarse grid 
% (note that matlab doesn't interpolate intelligently)
xx=-5:0.25:5; yy=-5:0.25:5;
xtest=[];
for i=1:length(xx)
    for j=1:length(yy)
        xtest=[xtest [xx(i) yy(j)]'];
    end
end

% yet another dataset for training and testing
N=[20,40];D=2;sz=[N(1),1];
tmp=randn(sz)+10;
angle=rand(sz)*2*pi-pi;
sumN=sum(sum(N));
X=[tmp.*cos(angle), tmp.*sin(angle);randn([N(2),D])]+...
  repmat((randn([1,D])+15),[sumN,1]);
Y=[zeros([N(1),1]);ones([N(2),1])];
idx=randperm(sum(N));
X=X(idx,:);Y=Y(idx,:);
xtrain=X';ctrain=Y';
N=sum(N);
minX=min(X);maxX=max(X);
bound=(maxX-minX)*0.1;
range_row=linspace(minX(1)-bound(1),maxX(1)+bound(1),30);
range_col=linspace(minX(2)-bound(2),maxX(2)+bound(2),30);
xx=range_row;
yy=range_col;
[X_map_test,Y_map_test]=meshgrid(range_row,range_col);
sumN_test=size(X_map_test,1)*size(X_map_test,2);
X_test = [reshape(X_map_test,[sumN_test,1]), ...
          reshape(Y_map_test,[sumN_test,1])];
xtest=X_test';
% size(xtrain),size(ctrain),size(xtest)

[meantrain,meantest,logpygx]=GPclass(xtrain,ctrain,xtest);
% plot(meantrain,'o'); 
% hold on; 	
% plot(ctrain,'r+'); title('training class labels and fit'); 
% drawnow; 
% hold off;

figure; hold on
for n=1:N
    if ctrain(n)==0;
        plot(xtrain(1,n),xtrain(2,n),'ro','markersize',5);
    else
        plot(xtrain(1,n),xtrain(2,n),'go','markersize',5);
    end
    if meantrain(n)<0.5;
        plot(xtrain(1,n),xtrain(2,n),'ro','markersize',8);
    else
        plot(xtrain(1,n),xtrain(2,n),'go','markersize',8);
    end
end
[cs,hnd]=contour(xx,yy,reshape(meantest,length(xx),length(yy)));
% clabel(cs,hnd,'fontsize',10); axis([-5 5 -5 5]);



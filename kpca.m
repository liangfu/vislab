function params = kpca(X,D_low,test_and_disp)
% KPCA - Kernel Principal Component Analysis

%% generate test data set
if nargin==0,
  testid=1;
  N=[20,20];
  D=2;
  if testid==1,
    sz=[N(1),1];
    tmp=randn(sz)+10;
    angle=rand(sz)*2*pi-pi;
    sumN=sum(sum(N));
    X=[tmp.*cos(angle), tmp.*sin(angle);randn([N(2),D])]+...
      repmat((randn([1,D])+15),[sumN,1]);
  elseif testid==2,
    sigma = rand([length(N),D]).*15.0+5.0;
    X=nan([sum(N),D]);
    for ii=1:length(N)
      X(sum(N(1:ii-1))+1:sum(N(1:ii)),:) = ...
          randn(N(ii),D)+repmat(sigma(ii,:),[N(ii),1]);  
    end
  end
  X=X(randperm(size(X,1)),:); % shuffle
  D_low=2;
  test_and_disp=1;
end
if nargin<3,if nargin<2,D_low=20;end,test_and_disp=1;end

%% shape of input data
sumN=size(X,1);
D=size(X,2);

%% using the Gaussian kernel to construct the kernel K
% K = zeros([sumN,sumN]);
% for r = 1:sumN
% for c = 1:r
% K(r,c) = exp(-sum(sum(((X(r,:) - X(c,:)).^2)))./2);
% K(c,r) = K(r,c);
% end
% end
K=exp(-0.5*(sum(X.^2,2)*ones(1,sumN)+ones(sumN,1)*sum(X.^2,2)'-2*(X*X')));

%% centering in feature space
unit = ones([sumN,sumN])./sumN;
Kc = K - unit*K - K*unit + unit*K*unit;

%% obtaining the low dimensional projection
[eigvec eigval]=eig(Kc);
eigval=real(diag(eigval));
% normalize
% for i=1:sumN, eigvec(:,i) = eigvec(:,i)/(sqrt(eigval(i))); end 
% eigvec=eigvec./repmat(sqrt(eigval(:))',[sumN,1]);
eigvec=eigvec.*repmat((eigval(:).^(-.5))',[sumN,1]);
[~, idx] = sort(eigval,'descend');
eigvec = eigvec(:,idx);

%% projecting the data in lower dimensions
Vk=eigvec(:,1:min(D_low,D));
beta = K*Vk; proj=beta;

%% pre-image computation
% DIST=zeros(sumN,sumN);
% for k=1:size(X,2), [a,b]=meshgrid(X(:,k)); DIST=DIST+(a-b).^2; end
% DIST=DIST.^0.5;
% DIST=DIST+diag(ones(sumN,1)*inf);
% para=10*median(min(DIST));

% maxiter=200;
% preimages=nan([sumN,D]);
% for item=1:sumN
%   gamma=zeros([1,sumN]);
%   for ii=1:sumN,gamma(ii) = Vk(ii,:)*beta(item,:)';end
%   z=mean(X)';
%   for iter=1:maxiter
%     z_old = z;
%     xx=bsxfun(@minus,X',z);
%     xx=xx.^2;
%     xx=-sum(xx)/(2*para.^2);
%     xx=exp(xx);
%     xx=xx.^gamma;
    
%     z=xx*X/sum(xx);
%     z=z';
%     if norm(z_old-z)/norm(z)<1e-5,break;end
%   end
%   preimages(item,:)=z;
% end

%% construct output parameters
params.proj=beta;
params.K=K;
params.eigvec=eigvec(:,1:min(D_low,D));
params.sumN=sumN;
params.X=X;

% test and display
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if test_and_disp

%% options
planar_disp=1;
nlevels = 20;

%% generate grid data for testing
minX=min(X);
maxX=max(X);
bound=(maxX-minX)*0.1;
range_row=linspace(minX(1)-bound(1),maxX(1)+bound(1),20);
range_col=linspace(minX(2)-bound(2),maxX(2)+bound(2),20);
[X_map_test,Y_map_test]=meshgrid(range_row,range_col);
sumN_test=size(X_map_test,1)*size(X_map_test,2);
X_test = [reshape(X_map_test,[sumN_test,1]), ...
          reshape(Y_map_test,[sumN_test,1]),zeros([sumN_test,size(X,2)-2])];
%% start testing 
% K_test = zeros([sumN_test,sumN]);
% for r = 1:sumN_test
% for c = 1:sumN
%   K_test(r,c) = exp(-0.5*sum(sum(((X_test(r,:) - X(c,:)).^2))));
% end
% end
K_test=exp(-0.5*(sum(X_test.^2,2)*ones(1,sumN)+...
                 ones(sumN_test,1)*sum(X.^2,2)'-2*(X_test*X')));

unit_test = ones([sumN_test,sumN])./sumN;
K_test_center=K_test-unit_test*K-K_test*unit+unit_test*K*unit;
proj_test=K_test_center*eigvec(:,1:min(D_low,D));

%% display
clf;
if 0,shading('faceted');else,shading('interp');end
colormap('gray'),

disptest=1;
%% show each component
for ii=1:2
lower=find(proj(:,ii)<=0);
upper=find(proj(:,ii)>0);
subplot(1,3,ii),
hold on,
if disptest
axis([minX(1)-bound(1),maxX(1)+bound(1),minX(2)-bound(2),maxX(2)+bound(2)]),
if planar_disp
  pcolor(range_row,range_col,reshape(proj_test(:,ii),size(X_map_test)))
else
  surf(range_row,range_col,reshape(proj_test(:,ii),size(X_map_test))),view(3)
end
end
% contour(range_row,range_col,...
%         reshape(proj_test(:,2),size(X_map_test)),nlevels,'b')
plot(X(lower,1),X(lower,2), 'r.',X(upper,1),X(upper,2), 'g.'),
title([num2str(ii) '-th component'])
hold off;
end

%% projection in feature space
lower=find(proj(:,1)<=0);
upper=find(proj(:,1)>0);
subplot(133),
plot(proj(lower,1),proj(lower,2),'r.',proj(upper,1),proj(upper,2),'b.'),
title('projection');

end % test_and_disp


function kpca2

N=[80,20];
D=2;
sz=[N(1),1];
tmp=randn(sz)+10;
angle=rand(sz)*2*pi-pi;
sumN=sum(sum(N));
X=[tmp.*cos(angle), tmp.*sin(angle);randn([N(2),D])]+...
  repmat((randn([1,2])+15),[sumN,1]);

%% using the Gaussian kernel to construct the kernel K
K = zeros([sumN,sumN]);
for r = 1:sumN
for c = 1:r
K(r,c) = exp(-sum(sum(((X(r,:) - X(c,:)).^2)))./2);
K(c,r) = K(r,c);
end
end

%% centering in feature space
unit = ones([sumN,sumN])./sumN;
K_center = K - unit*K - K*unit + unit*K*unit;

%% obtaining the low dimensional projection
[eigvec eigval]=eig(K_center);
eigval=real(diag(eigval));
% normalize
for i=1:sumN, eigvec(:,i) = eigvec(:,i)/(sqrt(eigval(i))); end 
[~, idx] = sort(eigval,'descend');
eigvec = eigvec(:,idx);

%% projecting the data in lower dimensions
proj = K_center*eigvec(:,1:D);

%% generate grid data for testing
x_range=0:2:30;
[X_map_test,Y_map_test]=meshgrid(x_range,x_range);
sumN_test=size(X_map_test,1)*size(X_map_test,2);
X_test = [reshape(X_map_test,[sumN_test,1]), ...
          reshape(Y_map_test,[sumN_test,1])];

%% start testing 
K_test = zeros([sumN_test,sumN]);
for r = 1:sumN_test
for c = 1:sumN
K_test(r,c) = exp(-sum(sum(((X_test(r,:) - X(c,:)).^2)))./2);
end
end
unit_test = ones([sumN_test,sumN])./sumN;
K_test_center=K_test-unit_test*K-K_test*unit+unit_test*K*unit;
proj_test=K_test_center*eigvec(:,1:D);

%% display
clf;
shading('faceted'),
colormap('gray'),

%% first component
subplot(131),
hold on,
pcolor(x_range,x_range,reshape(proj_test(:,1),size(X_map_test))),
contour(x_range,x_range,reshape(proj_test(:,1), ...
                                        size(X_map_test)),12,'b')
plot(X(1:N(1),1),X(1:N(1),2), 'b.',...
     X(N(1)+1:N(1)+N(2),1),...
     X(N(1)+1:N(1)+N(2),2), 'r.'),...
    title('original data'),
hold off;

%% second component
subplot(132),
hold on,
pcolor(x_range,x_range,reshape(proj_test(:,2),size(X_map_test))),
contour(x_range,x_range,reshape(proj_test(:,2), ...
                                       size(X_map_test)),12,'b')
plot(X(1:N(1),1),X(1:N(1),2), 'b.',...
     X(N(1)+1:N(1)+N(2),1),...
     X(N(1)+1:N(1)+N(2),2), 'r.'),...
    title('original data'),
hold off;

%% projection in feature space
subplot(133),plot(proj(1:N(1),1),proj(1:N(1),2),'b.', ...
                  proj(N(1)+1:N(1)+N(2),1),...
                  proj(N(1)+1:N(1)+N(2),2),'r.'), ...
    title('projection');


function chanvese
sigma=3;
[I,C,phi]=testimg();
k=fspecial('gaussian',[5,5],3);I=imfilter(I,k); % gaussian smooth 

for iter=1:1000
% compute dirac term
diracPhi=Dirac(phi,sigma);

% compute curvature
[dx,dy]=gradient(phi);
s=sqrt(dx.^2+dy.^2+1e-5);
nx=dx./s;ny=dy./s;

% compute mu inside/outside
Hv=Heaviside(phi,sigma);
mu_ex=sum(sum(I.*Hv))/sum(Hv(:));
mu_in=sum(sum(I.*(1-Hv)))/sum(1-Hv(:));

% minimize Chan-Vese energy 
dphi=diracPhi.*(0.02*div(nx,ny)-(I-mu_ex).^2+(I-mu_in).^2);
dt = .48/(max(abs(dphi(:)))+1e-5);

phi=phi+dt.*dphi;
phi=sussman(phi,0.5); % reinitialize

% termination criteria
mask_new=phi<0;
if ~exist('mask_old'),mask_old=zeros(size(mask_new));end
mask_diff=abs(mask_new-mask_old);
if sum(mask_diff(:))<2,disp(sprintf('converged within %d iter',iter));break;end
mask_old=mask_new;

% display
if mod(iter,20)==0,
boundaries=bwboundaries(phi<0);if length(boundaries)<1,break;end;
subplot(121),imshow(I,[]);hold on;
for ii=1:length(boundaries),C=boundaries{ii};plot(C(:,2),C(:,1),'r-','LineWidth',2);end;hold off;
subplot(122),imagesc(phi),title(sprintf('%f,',[mu_in,mu_ex]));
drawnow;pause(.1)
end

end % end of iter loop

function f = div(nx,ny)
[nxx,~]=gradient(nx);  
[~,nyy]=gradient(ny);
f=nxx+nyy;

function f = Dirac(x, sigma)
f=(0.5/sigma)*(1+cos(pi*x/sigma));
b = (x<=sigma) & (x>=-sigma);
f = f.*b;

function H = Heaviside(phi, sigma)
idx0 = find((phi<=sigma) & (phi>=-sigma));
idx1 = find(phi<-sigma);
idx2 = find(phi>sigma);
H = zeros(size(phi,1),size(phi,2));
H(idx0) = ((0.5/sigma)*phi(idx0) + 0.5/pi*sin(pi*phi(idx0)/sigma))+0.5; 
H(idx1)=1e-5;                                  
H(idx2)=1-1e-5;                                 

function D = sussman(D, dt)
% forward/backward differences
a = D - shiftR(D); % backward
b = shiftL(D) - D; % forward
c = D - shiftD(D); % backward
d = shiftU(D) - D; % forward
a_p = a;  a_n = a; % a+ and a-
b_p = b;  b_n = b;
c_p = c;  c_n = c;
d_p = d;  d_n = d;
a_p(a < 0) = 0; a_n(a > 0) = 0;
b_p(b < 0) = 0; b_n(b > 0) = 0;
c_p(c < 0) = 0; c_n(c > 0) = 0;
d_p(d < 0) = 0; d_n(d > 0) = 0;
dD = zeros(size(D));
D_neg_ind = find(D < 0);
D_pos_ind = find(D > 0);
dD(D_pos_ind) = sqrt(max(a_p(D_pos_ind).^2, b_n(D_pos_ind).^2) ...
                     + max(c_p(D_pos_ind).^2, d_n(D_pos_ind).^2)) - 1;
dD(D_neg_ind) = sqrt(max(a_n(D_neg_ind).^2, b_p(D_neg_ind).^2) ...
                     + max(c_n(D_neg_ind).^2, d_p(D_neg_ind).^2)) - 1;
D = D - dt .* sussman_sign(D) .* dD;
function shift = shiftD(M)
shift = shiftR(M')';
function shift = shiftL(M)
shift = [ M(:,2:size(M,2)) M(:,size(M,2)) ];
function shift = shiftR(M)
shift = [ M(:,1) M(:,1:size(M,2)-1) ];
function shift = shiftU(M)
shift = shiftL(M')';
function S = sussman_sign(D)
S = D ./ sqrt(D.^2 + 1);   

function [img,c,phi]=testimg
% generate fixed sized test image
bw=zeros([256,256]);
nr=size(bw,1);
nc=size(bw,2);
noisy=int32(nc*nr*.3);
bw(int32(nr/4):int32(nr/4*3),int32(nc/4):int32(nc/4*3))=1;
x=randi([1,255],[noisy 2]);y=x(:,1)*256+x(:,2);bw(y)=1;
x=randi([1,255],[noisy 2]);y=x(:,1)*256+x(:,2);bw(y)=0;
img=double(bw);
% generate initial curve
angle=1e-5:.2:pi*2.;
tmp=randi([38,48],[1,max(size(angle))]);
c=int32([tmp.*cos(angle)+80;tmp.*sin(angle)*1.8+120]');
% display them, optional
%imshow(img,[]);hold on,plot(c(:,1),c(:,2),'r-','LineWidth',2),hold off;
c=double(c);
mask=zeros(size(bw));
mask=roipoly(mask,c(:,1),c(:,2));
phi=bwdist(mask)-bwdist(mask==0);


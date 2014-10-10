function scwad
% SCWAD - Shadow Compensation using Weighted Average Difference

clear all;clf;
path='../data/illumination/';
a=dir([path '*.pgm']);
N=length(a);
for ii=1:N
im = double(imread([path a(ii).name]))./256.;
idx=find(im<=0);im(idx)=1e-5;

% gamma correction
re=im.^.2;

% DoG filtering
sigma0=1.0;sigma1=sigma0*2.;
kernel0=fspecial('gaussian',sigma0*3.,sigma0);
kernel1=fspecial('gaussian',sigma1*3.,sigma1);
g0=imfilter(re,kernel0);
g1=imfilter(re,kernel1);
re=g0-g1;

re=re./(mean(abs(re(:)).^.1).^10.);
% re=re./(mean(min(abs(re(:))).^.1).^10.);
re=tanh(re./10.)*10.;

% minval=min(min(re));maxval=max(max(re));% [minval,maxval]
% re=re-minval;re=re./(maxval-minval);
% re=re./std(re(:));
% mu=mean(re(:)); 
% re=(re-mu+.5);
% idx=find(re<-1.5);re(idx)=-1.5;
% idx=find(re>7.5);re(idx)=7.5;

% minval=min(min(re));maxval=max(max(re));[minval,maxval]
% re=(re-minval)./(maxval-minval)+1e-5;
% re=re.^.8;

% re=histeq(re);

subplot(N,2,1+(ii-1)*2),imshow(im,[]);
subplot(N,2,2+(ii-1)*2),imshow(re,[]);
end







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% beta=[0.3,2.,3.];

% A=(-5:.1:5)';N=length(X);
% y=(A.^2)*beta(1)+A*beta(2)+randn([N,1])*.2+beta(3);
% idx=randperm(N);
% A=A(idx);y=y(idx);
% invS=cov(y)^-1;
% x=inv(A'*invS*A)*A'*invS*y;
% x=inv(A'*A)*A'*y;

% clf;
% plot(A,y,'.');

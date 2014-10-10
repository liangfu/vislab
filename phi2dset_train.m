function phi2dset_train

clear all;
fp=fopen('../data/phi2dset.txt','r');
szstr=fgets(fp);
sz=sscanf(szstr,'%d %d %d\n');
N=sz(1);nr=sz(2);nc=sz(3);[N,nr,nc]
dat=fread(fp,N*nr*nc,'single');
% dat = reshape(dat,[],N);
dat = reshape(dat,[],N);
% imshow(idct2(reshape(dat(:,1),16,16)'),[]);

% append training data set
if 1,
wdir='../data/close/';
files = dir([wdir '*.bmp']);
for ii=1+N:length(files)+N
  bw=imread([wdir files(ii-N).name])>0;
  im_rst = uint8(imnormalize(bw,size(bw,1)*size(bw,2)*0.25));
  phi = double(initsdf(im_rst));
  phi_dct=dct2(phi);
  phi_dct=phi_dct(1:16,1:16);
  dat(:,ii)=reshape(idct2(phi_dct)',[],1);
end
N=size(dat,2);
end

newsz = [16,16];

% X=dat';
for ii=1:N
tmp=reshape(dat(:,ii),nc,nr)';
% tmp=dct2(tmp); tmp=tmp(1:newsz(1),1:newsz(2));
% imgsz=[40,30];
% tmp=idct2([[tmp,zeros(newsz(1),imgsz(2)-newsz(2))];...
%            zeros(imgsz(1)-newsz(1),imgsz(2))]);
% imshow(tmp,[]);
% tmp=idct2(tmp);
tmp=tmp';
X(ii,:)=reshape(tmp,1,[]);
end
% imshow(reshape(mean(X),16,16)',[]);
% X(1,1:5)

% nframes=8;
% for ii=1:N-(nframes-1),
% stip(ii,:)=reshape(X(ii:ii+(nframes-1),:),1,[]);
% end

Y=zeros([N,1]);
Y(329:522)=1;Y(562:589)=1;Y(656:579)=1;Y(796:843)=1;Y(851:906)=1;
Y(1:5)=-1;Y(97:104)=-1;Y(491:501)=-1;Y(565:571)=-1;Y(845:850)=-1;

trnidx=find(Y>=0);
X=X(trnidx,:);Y=Y(trnidx);

% [W,PRIOR]=kfd(X(100:end-20,:),Y(100:end-20,:));
[W,PRIOR]=lda(X(end-200:end,:),Y(end-200:end));
sprintf('%ff,',W')
sprintf('%ff,',PRIOR')

opt = kpca(X,int32(size(X,2)),0); proj=opt.proj;
% [~,~,~,proj]=pca(X);
% W=kfd(real(proj),Y);
% W=kfd(X,Y);
% W=lda(X,Y);sprintf('%.5f,',W)
W=lda(real(proj),Y);
% W=lda_new(real(opt.proj),Y);

startp=1+49*4;
clf
for ii=startp:min(N,startp+48)
subplot(7,7,ii-startp+1),
imshow(reshape(X(ii,:),nc,nr)'>0),
title(num2str(ii))%,drawnow
end

% size(dat)
% mu(:,:)=mean(dat,3); 
% mu=mu';
% imagesc(mu);
% fclose(fp);


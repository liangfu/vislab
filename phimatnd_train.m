function phimatnd_train

clear all;
fp=fopen('../data/phimatnd.txt','r');
szstr=fgets(fp);
sz=sscanf(szstr,'%d %d %d\n');
N=sz(1);nr=sz(2);nc=sz(3);[N,nr,nc]
dat = reshape(fread(fp,N*nr*nc,'single'),[],N);

newsz = [16,16];

% X=dat';
for ii=1:N
tmp=reshape(dat(:,ii),nc,nr)';
tmp=dct2(tmp); tmp=tmp(1:newsz(1),1:newsz(2));
% imgsz=[40,30];
% tmp=idct2([[tmp,zeros(newsz(1),imgsz(2)-newsz(2))];...
%            zeros(imgsz(1)-newsz(1),imgsz(2))]);
% imshow(tmp,[]);
tmp=idct2(tmp);tmp=tmp';
X(ii,:)=reshape(tmp,1,[]);
end

% nframes=8;
% for ii=1:N-(nframes-1),
% stip(ii,:)=reshape(X(ii:ii+(nframes-1),:),1,[]);
% end

Y=zeros([N,1]);
Y(329:522)=1;Y(560:587)=1;Y(657:577)=1;Y(796:843)=1;

W=lda_new(X(100:end-20,:),Y(100:end-20,:));
sprintf('%.2f,',W)

opt = kpca(X,int32(size(X,2)),0); proj=opt.proj;
% [~,~,~,proj]=pca(X);
% W=kfd(real(proj),Y);
% W=kfd(X,Y);
% W=lda_new(X,Y);
W=lda(real(proj),Y);
% W=lda_new(real(opt.proj),Y);


% startp=1+100*8;
% clf
% for ii=startp:min(N,startp+99)
% subplot(10,10,ii-startp+1),
% imshow(reshape(dat(:,ii),nc,nr)'>0),
% title(num2str(ii)),drawnow
% end

% size(dat)
% mu(:,:)=mean(dat,3); 
% mu=mu';
% imagesc(mu);
% fclose(fp);


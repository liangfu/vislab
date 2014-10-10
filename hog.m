function [H,G]=hog(I,csize,bound)
% HOG function extracts histogram of oriented gradients
% Example :
%    I=rgb2gray(imread('../data/lena.jpg'));
%    cellsize=6;
%    [H,G]=hog(I,cellsize)
% Or,
%    clf;[~,G]=hog(gabor2,10,1);
%    subplot(121),imagesc(gabor2),
%    subplot(122),imagesc(G);

if nargin<3,bound=1;
if nargin<2,csize=6;
if nargin<1,I=gabor2;
end;end;end

nbins=9;

if size(I,3)~=1,I=rgb2gray(I);end
I=histeq(I);
I=double(I);

hx=[-1,0,1];
hy=hx';
dx=imfilter(I,hx);
dy=imfilter(I,hy);

mag=sqrt(dx.^2+dy.^2);
ang=atan2(dy,dx);

% idx=find(ang<0);
% ang(idx)=ang(idx)+pi;
% ang=int32(ang/pi*8+1);

% clf;
% subplot(141),imagesc(mag);
% subplot(142),imagesc(ang/pi*180),title('as is');
% loweridx=find(ang<0);
% ang(loweridx)=-1;
% subplot(143),imagesc(ang/pi*180),title('set to -1');
% loweridx=find(ang<0);
% ang(loweridx)=ang(loweridx)+pi;
% subplot(144),imagesc(ang/pi*180),title('+pi');
% drawnow;

cont=0;
for bm=1:(size(I,1)-2*bound)/csize
for bn=1:(size(I,2)-2*bound)/csize
bmpos=(bm-1)*csize+bound+1:csize+(bm-1)*csize+bound;
bnpos=(bn-1)*csize+bound+1:csize+(bn-1)*csize+bound;
subang = ang(bmpos,bnpos);
submag = mag(bmpos,bnpos);
subang = reshape(subang,[],1);
submag = reshape(submag,[],1);
bin=1;Hcurr=zeros([nbins,1]);
idx=find(subang<0);subang(idx)=subang(idx)+pi;subang=int32(subang/pi*8+1);
for hval=1:9
  idx=find(subang==hval);
  Hcurr(bin)=sum(submag(idx));
  bin=bin+1;
end
Hcurr=Hcurr./(norm(Hcurr)+0.01);
H(cont*nbins+1:cont*nbins+nbins)=Hcurr;
cont=cont+1;
end
end

G=zeros(size(I));
[X,Y]=meshgrid(-(csize-1)/2:(csize-1)/2,-(csize-1)/2:(csize-1)/2);
cont=0;
for bm=1:(size(I,1)-2*bound)/csize
for bn=1:(size(I,2)-2*bound)/csize
bmpos=(bm-1)*csize+bound+1:(bm-1)*csize+csize+bound;
bnpos=(bn-1)*csize+bound+1:(bn-1)*csize+csize+bound;
hvals=H(cont*nbins+1:cont*nbins+nbins);
for bin=1:nbins
  G(bmpos,bnpos)=G(bmpos,bnpos)+...
      gauss2d(X,Y,pi*(bin-1)/9,[csize*2.,csize/20]).*hvals(bin);
      % gauss2d(X,Y,pi*(bin-1)/9,[csize/20,csize*2.]).*hvals(bin);
end
cont=cont+1;
end
end



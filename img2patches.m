function D=img2patches(im,winsize)

nr=floor(size(im,1)/winsize);
nc=floor(size(im,2)/winsize);
% [X,Y]=meshgrid(1:winsize:size(im,2)-winsize+1,...
%                1:winsize:size(im,1)-winsize+1);
[X,Y]=meshgrid(1:winsize:nr*winsize,...
               1:winsize:nc*winsize);
imgloc=[reshape(X,[],1) reshape(Y,[],1)];
% imgloc
for ii=1:size(imgloc,1)
sample=im(imgloc(ii,1):imgloc(ii,1)+winsize-1,...
          imgloc(ii,2):imgloc(ii,2)+winsize-1);
D(ii,:)=reshape(sample,1,[]);
end
function im2=patches2img(samples,winsize,imsize)

im2=zeros(imsize);
nr=floor(imsize(1)/winsize);
nc=floor(imsize(2)/winsize);
% [X,Y]=meshgrid(1:winsize:size(im,2)-winsize+1,...
%                1:winsize:size(im,1)-winsize+1);
[X,Y]=meshgrid(1:winsize:nr*winsize,...
               1:winsize:nc*winsize);
imgloc=[reshape(X,[],1) reshape(Y,[],1)];

for ii=1:size(samples,1)
im2(imgloc(ii,1):imgloc(ii,1)+winsize-1,...
    imgloc(ii,2):imgloc(ii,2)+winsize-1)=reshape(samples(ii,:),winsize,[]);
end

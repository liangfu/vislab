function fval=iconicfeature(img)
% ICONICFEATURE extract 27 iconic feature of an image patch

% clear all; clf;
img=imresize(img,32.0/size(img,1));
    
ii=1;fval=[];
for sigma=[1,2,4]
sz=sigma*8;
[X,Y]=meshgrid(-(sz-1)/2.:(sz-1)/2.);
G0=exp(-(X.^2+Y.^2)./(2*sigma^2));
G1=imfilter(G0,[-1,0,1]);
G2=imfilter(G0,[-1;0;1]);
G3=imfilter(G1,[-1,0,1]);
G4=imfilter(G1,[-1,0,0;0,0,0;0,0,1]);
G5=imfilter(G1,[0,0,1;0,0,0;-1,0,0]);
G6=imfilter(G3,[-1,0,1]);
G7=imfilter(G4,[-1,0,0;0,0,0;0,0,1]);
G8=imfilter(G5,[0,0,-1;0,0,0;1,0,0]);
G9=imfilter(G2,[-1;0;1]);
G9=imfilter(G9,[-1;0;1]);

fval(ii)=sum(sum(imfilter(img,G1,'replicate'))); ii=ii+1;
fval(ii)=sum(sum(imfilter(img,G2,'replicate'))); ii=ii+1;
fval(ii)=sum(sum(imfilter(img,G3,'replicate'))); ii=ii+1;
fval(ii)=sum(sum(imfilter(img,G4,'replicate'))); ii=ii+1;
fval(ii)=sum(sum(imfilter(img,G5,'replicate'))); ii=ii+1;
fval(ii)=sum(sum(imfilter(img,G6,'replicate'))); ii=ii+1;
fval(ii)=sum(sum(imfilter(img,G7,'replicate'))); ii=ii+1;
fval(ii)=sum(sum(imfilter(img,G8,'replicate'))); ii=ii+1;
fval(ii)=sum(sum(imfilter(img,G9,'replicate'))); ii=ii+1;
end

fval=fval./(32.0*32.0);

% subplot(2,5,1),imshow(G1,[]);
% subplot(2,5,2),imshow(G2,[]);
% subplot(2,5,3),imshow(G3,[]);
% subplot(2,5,4),imshow(G4,[]);
% subplot(2,5,5),imshow(G5,[]);
% subplot(2,5,6),imshow(G0,[]);
% subplot(2,5,7),imshow(G6,[]);
% subplot(2,5,8),imshow(G7,[]);
% subplot(2,5,9),imshow(G8,[]);
% subplot(2,5,10),imshow(G9,[]);

% subplot(2,5,1),imshow(imfilter(img,G1,'replicate'),[]);
% subplot(2,5,2),imshow(imfilter(img,G2,'replicate'),[]);
% subplot(2,5,3),imshow(imfilter(img,G3,'replicate'),[]);
% subplot(2,5,4),imshow(imfilter(img,G4,'replicate'),[]);
% subplot(2,5,5),imshow(imfilter(img,G5,'replicate'),[]);
% subplot(2,5,6),imshow(imfilter(img,G0,'replicate'),[]);
% subplot(2,5,7),imshow(imfilter(img,G6,'replicate'),[]);
% subplot(2,5,8),imshow(imfilter(img,G7,'replicate'),[]);
% subplot(2,5,9),imshow(imfilter(img,G8,'replicate'),[]);
% subplot(2,5,10),imshow(imfilter(img,G9,'replicate'),[]);

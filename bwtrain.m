function bwtrain

nr=50;nc=50;
rows=16;cols=24;

wdir='../data/palmdata/palmdata1/';
files=dir([wdir '*.bmp']);

X=zeros([length(files),rows*cols]);
for ii=1:100
bw=imread([wdir files(ii).name])>0;

% extract features
im_rst = uint8(imnormalize(bw,size(bw,1)*size(bw,2)*0.25));
phi = double(initsdf(im_rst));
phi_dct=dct2(phi);
phi_dct=phi_dct(1:rows,1:cols);
X(ii,:)=reshape(idct2(phi_dct)',1,[]);
end

N=inf;
startp=1+49*1;
clf
for ii=startp:min(N,startp+48)
subplot(7,7,ii-startp+1);
imshow(reshape(X(ii,:),cols,rows)'>0,[]);
end

% for ii=1:49
%   subplot(7,7,mod(ii,49));
%   imshow(reshape(X(ii,:),rows,cols)'>0,[]);
% end

Y=zeros([size(X,1),1]);
Y(1:7)=1;

[W,PRIOR]=lda(X,Y);

function x=extractfeature(bw,objsz)
im_rst = uint8(imnormalize(bw,objsz));
phi = double(initsdf(im_rst));
phi_dct=dct2(phi);
phi_dct=phi_dct(1:rows,1:cols);
% imshow(idct2(phi_dct),[]);
x=reshape(phi_dct',1,[]);

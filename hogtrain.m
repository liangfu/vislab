function hogtrain
% HOGTRAIN train images with hog feature

clear all;
wdir='../data/hogtrain/RotateSeq/';
files=dir([wdir '*.pgm']);
bound=1;
csize=6;
% wdir='../data/hogtrain/pedestrians128x64/';
% files=dir([wdir '*.ppm']);

% Tsz=32;

for m=1:10:length(files)
% I=imresize(imread([wdir files(m).name]),[26,26],'cubic');
I=imread([wdir files(m).name]);
I=I(3:end-2,:,:);
%I=imresize(I,[Tsz,Tsz],'cubic');
I=I(5:end-4,5:end-4);
% [H(m,:),G]=hog(I,8,4);
[Hpos(floor(m/10)+1,:),G]=hog(I,csize,bound);break;
end
Tsz=size(I,1)

clf;
subplot(121),imshow(I,[]);
subplot(122),imshow(G,[]);

wdir='../data/hogtrain/scene/';
files=dir([wdir '*.png']);
for m=1:size(Hpos,1)*2
I=imread([wdir files(randi(length(files))).name]);
Isz=size(I);Isz=Isz(1:2);
loc=randi([bound,min(Isz)-bound-Tsz],[1,2]);
I=I(loc(1):loc(1)+Tsz-1,loc(2):loc(2)+Tsz-1,:);
[Hneg(m,:),G]=hog(I,csize,bound);
end

X=[Hpos;Hneg];
Y=zeros([size(Hpos,1)+size(Hneg,1),1]);
Y(1+size(Hpos,1):end)=1;
[W,PRIOR]=lda(X,Y,min(size(X)));

I=imread('../data/hogtrain/raw/00180.png');
P=hogdetect(I,W,PRIOR,csize,bound);

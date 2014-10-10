function partstruct
% PARTSTRUCT builds pictorial structure model
% 
% See also
%    getsubrect
% 
% Copyright (c) Liangfu Chen (2012-2013)


partsName={'leye','reye','nose','lmou','cmou','rmou'};

clear all;
% wdir='../data/face/';
% ddir='images/';
% gtruthfn='ground_truth_a.txt';
wdir='../dataset/yalefaces/';
ddir='';
gtruthfn='annotation.txt';

fp = fopen([wdir gtruthfn],'r');
gtruth=[];
I=1;N=0;
while 1
  gtruth(I).name = fscanf(fp, '%s',1); 
  if length(gtruth(I).name)==0 || gtruth(I).name(1)=='-', N=I-1;break; end
  % gtruth(I).data = fscanf(fp, '%f',12);
  gtruth(I).data = fscanf(fp, '%f',10);
  gtruth(I).data = reshape(gtruth(I).data,2,[])';
  I=I+1;
end
fclose(fp);

sz=32;
for I=1:N
data=gtruth(I).data;
scale=sqrt((data(1,1)-data(2,1)).^2+(data(1,2)-data(2,2)).^2)/(sz*1.5);
img = imread([wdir ddir gtruth(I).name]);
leye(I,:,:)=imresize(getsubrect(img,data(1,:),[scale*sz,scale*sz]),[sz,sz]);
reye(I,:,:)=imresize(getsubrect(img,data(2,:),[scale*sz,scale*sz]),[sz,sz]);
nose(I,:,:)=imresize(getsubrect(img,data(3,:),[scale*sz,scale*sz]),[sz,sz]);
lmou(I,:,:)=imresize(getsubrect(img,data(4,:),[scale*sz,scale*sz]),[sz,sz]);
%mmou(I,:,:)=imresize(getsubrect(img,data(5,:),[scale*sz,scale*sz]),[sz,sz]);
rmou(I,:,:)=imresize(getsubrect(img,data(5,:),[scale*sz,scale*sz]),[sz,sz]);
end

% subplot(331),imshow(imresize(leye,2),[]); colormap gray;
% subplot(333),imshow(imresize(reye,2),[]); colormap gray;
% subplot(335),imshow(imresize(nose,2),[]); colormap gray;
% subplot(337),imshow(imresize(lmou,2),[]); colormap gray;
% subplot(338),imshow(imresize(mmou,2),[]); colormap gray;
% subplot(339),imshow(imresize(rmou,2),[]); colormap gray;

clf;
subplot(331),imshow(imresize(reshape(mean(leye,1),sz,sz),2),[]); 
subplot(333),imshow(imresize(reshape(mean(reye,1),sz,sz),2),[]); 
subplot(335),imshow(imresize(reshape(mean(nose,1),sz,sz),2),[]); 
subplot(337),imshow(imresize(reshape(mean(lmou,1),sz,sz),2),[]); 
% subplot(338),imshow(imresize(reshape(mean(mmou,1),sz,sz),2),[]); 
subplot(339),imshow(imresize(reshape(mean(rmou,1),sz,sz),2),[]); 

for I=1:N
fval(1,I,:)=iconicfeature(reshape(leye(I,:,:),sz,sz));
fval(2,I,:)=iconicfeature(reshape(reye(I,:,:),sz,sz));
fval(3,I,:)=iconicfeature(reshape(nose(I,:,:),sz,sz));
fval(4,I,:)=iconicfeature(reshape(lmou(I,:,:),sz,sz));
fval(5,I,:)=iconicfeature(reshape(rmou(I,:,:),sz,sz));
end

for I=1:5
parts(I).mean=mean(reshape(fval(I,:,:),[],27),1);
parts(I).covar=cov(reshape(fval(I,:,:),[],27),1);
end



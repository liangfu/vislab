function dictcollect

clear all;clf;
% wdir='../dataset/Train/neg/';
wdir='../dataset/upperbodyfrontal_dataset/images/';
% wdir='../dataset/CroppedYale/yaleB24/';
files=dir(wdir);
imgidx=randi([3,length(files)],[11000,1]);
for ii=1:length(imgidx)
im=imread([wdir files(imgidx(ii)).name]);
if size(im,3)==3,im=rgb2gray(im);end
imgloc=[randi(size(im,1)-7) randi(size(im,2)-7)];
sample=double(im(imgloc(1):imgloc(1)+7,imgloc(2):imgloc(2)+7));
D(ii,:)=reshape(sample,1,[]);
% if ii==500,break;end
end

%% sort sample dictionary element by variance
%% display the 500 dictionary elements
% [~,idx]=sort(std(D'));
% D=D(idx,:);
% nr=10;nc=50;
% disp=zeros([nr*8,nc*8]);
% for ii=1:nr
% for jj=1:nc
% disp((ii-1)*8+1:(ii-1)*8+8,(jj-1)*8+1:(jj-1)*8+8)=...
%     reshape(D((ii-1)*nc+jj,:),8,[]);
% end
% end
% imshow(disp,[]);

%% collect dictionary directly from face images
% load ../dataset/olivettifaces
% imgidx=randi(size(faces,2),[11000,1]);
% imgloc=randi(64-7,[length(imgidx) 2]);
% D=[];
% for ii=1:length(imgidx)
% im=reshape(faces(:,imgidx(ii)),64,[]);
% sample=im(imgloc(ii,1):imgloc(ii,1)+7,imgloc(ii,2):imgloc(ii,2)+7);
% D(ii,:)=reshape(sample,1,[]);
% % if ii==500,break;end
% end

size(D)
% assert(size(D,1)==11000);
assert(size(D,2)==64);

fp=fopen('../data/dictinit11000x64.bin','w');
fwrite(fp,D');
fclose(fp);

%% sort sample dictionary element by variance
%% display the 500 dictionary elements
[~,idx]=sort(std(D'));
D=D(idx,:);
nr=10;nc=50;
disp=zeros([nr*8,nc*8]);
for ii=1:nr
for jj=1:nc
disp((ii-1)*8+1:(ii-1)*8+8,(jj-1)*8+1:(jj-1)*8+8)=...
    reshape(D((ii-1)*nc+jj,:),8,[]);
end
end
imshow(disp,[]);


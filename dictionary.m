function dictionary

load ../dataset/olivettifaces
% size(faces)

imgidx=randi(size(faces,2),[11000,1]);
imgloc=randi(64-7,[length(imgidx) 2]);

D=[];

for ii=1:length(imgidx)
im=reshape(faces(:,imgidx(ii)),64,[]);
sample=im(imgloc(ii,1):imgloc(ii,1)+7,imgloc(ii,2):imgloc(ii,2)+7);
D(ii,:)=reshape(sample,1,[]);
% D(ii,:)=D(ii,:)-min(D(ii,:))+1;
% imshow(im,[]);drawnow;
if ii==500,break;end
end

size(D)
% assert(size(D,1)==11000);
assert(size(D,2)==64);

% fp=fopen('../data/dictinit11000x64.bin','w');
% fwrite(fp,D');
% fclose(fp);

%% sort sample dictionary element by variance
[~,idx]=sort(std(D'));
D=D(idx,:);

%% display the 500 dictionary elements
nr=10;nc=50;
disp=zeros([nr*8,nc*8]);
for ii=1:nr
for jj=1:nc
disp((ii-1)*8+1:(ii-1)*8+8,(jj-1)*8+1:(jj-1)*8+8)=...
    reshape(D((ii-1)*nc+jj,:),8,[]);
end
end

imshow(disp,[]);
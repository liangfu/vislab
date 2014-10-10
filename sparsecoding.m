function out=sparsecoding(im)

% clear all;clf;
if nargin<1,
% im=imread('../dataset/CroppedYale/yaleB01/yaleB01_P00A-005E-10.pgm');
im=imread('../dataset/olivettifaces/olivetti-001.png');
% im=imread('../dataset/CroppedYale/yaleB11/yaleB11_P00A-005E-10.pgm');
end
winsize=8;
alg=@OrthogonalMatchingPursuit;
% alg=@BasisPursuit;

%% load global dictionary vocabularies
% fp=fopen('../data/dictinit11000x64.bin','r');
% D=reshape(fread(fp,11000*(winsize^2),'uint8'),(winsize^2),[])';
% fclose(fp);

fp=fopen('../data/dict90x64.bin','r');
D=reshape(fread(fp,90*(winsize^2),'uint8'),(winsize^2),[])';fclose(fp);

%% get sub-windows & construct sparse representation for the signal
freq=zeros([size(D,1),1]);
samples=double(img2patches(im,winsize));% size(samples)
% size(samples,1),size(im)
for ii=1:size(samples,1)
out(ii)=alg(double(D'),samples(ii,:),8);
freq(out(ii).loc)=freq(out(ii).loc)+1;
end

%% record frequently used dictionary elements
% idx=find(freq>15);
% fp=fopen('../data/dict140x64.bin','w');fwrite(fp,D(idx,:)');fclose(fp);

%% print top dictionary elements
% [val,idx]=sort(freq,'descend');
% % val(1:10),idx(1:10)
% D=D(idx(1:100),:);

%% display dictionary elements
% %% sort sample dictionary element by variance
% % [~,idx]=sort(std(D')); D=D(idx,:);
% nr=10;nc=10;
% disp=zeros([nr*winsize,nc*winsize]);
% for ii=1:nr
% for jj=1:nc
% disp((ii-1)*winsize+1:(ii-1)*winsize+winsize,...
%      (jj-1)*winsize+1:(jj-1)*winsize+winsize)=...
%     reshape(D((ii-1)*nc+jj,:),winsize,[]);
% end
% end
% % subplot(141),
% imshow(disp,[]);

%% signal reconstruction
% for ii=1:size(samples,1)
% samples2(ii,1:winsize^2)=D(out(ii).loc,:)'*out(ii).val;
% end
% im2=patches2img(samples2,winsize,size(im));
% subplot(142),imshow(im,[]);
% subplot(143),imshow(im2,[]);
% subplot(144),imshow(sqrt((double(im)-im2).^2),[]);


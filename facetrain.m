function facetrain

clear all;clf;
Nper=16;
Nimg=1;
winsize=8;

%% load static dictionary
fp=fopen('../data/dict90x64.bin','r');
D=reshape(fread(fp,90*(winsize^2),'uint8'),(winsize^2),[])';fclose(fp);

%% pick random person and its random profile
perids=randperm(32);perids=perids(1:Nper)';
imgids=randi([1,5],[Nper,1]);%randi(10,[Nper,1]);

for ii=1:Nper
for jj=1:Nimg
id=(perids(ii)-1)*10+imgids(ii);
im=imread(['../dataset/olivettifaces/olivetti-' ...
           sprintf('%03d',id) '.png']);im=im(:,5:end-4);
% im=imread(['../dataset/CroppedYale/yaleB' sprintf('%02d',ii) ...
%            '/yaleB' sprintf('%02d',ii) '_P00A-005E-10.pgm']);   
out=sparsecoding(im);
aa(ii).id=id;aa(ii).im=im;aa(ii).out=out;
aa(ii).val=zeros(winsize^2,size(D,1));
for kk=1:size(out,2)
aa(ii).val(kk,out(kk).loc)=out(kk).val;
end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for testNimgs=1:20
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% select random individual
id=randi([331,400],1);
id=(perids(randi(Nper,1))-1)*10+randi([6,10],1);
testimg=imread(['../dataset/olivettifaces/olivetti-' ...
                sprintf('%03d',id) '.png']);testimg=testimg(:,5:end-4);
% testimg=imread(['../dataset/CroppedYale/yaleB' sprintf('%02d',id) ...
%                 '/yaleB' sprintf('%02d',id) '_P00A-005E-10.pgm']);   
tsamples=img2patches(testimg,winsize);

%% test image for face recognition
for ii=1:Nper
out=aa(ii).out;
for jj=1:length(out)
d=double(D(out(jj).loc,:))';
y=double(tsamples(jj,:))';
% x=d\y;%inv(d'*d+1e-5)*(d'*y); % least square
bpout=BasisPursuit(d,y',4);x=bpout.val;
bb(ii).out(jj).val=x;
end
% compute error
samples=[];
for jj=1:length(out)
% samples2(jj,1:winsize^2)=D(out(jj).loc,:)'*out(jj).val;
samples2(jj,1:winsize^2)=D(out(jj).loc,:)'*bb(ii).out(jj).val;
end
im2=patches2img(samples2,winsize,size(testimg));
bb(ii).im=im2;
% bb(ii).err=sqrt((im2-double(aa(ii).im)).^2);
bb(ii).err=sqrt((im2-double(testimg)).^2); %abs(im2-double(testimg));%
err(ii)=sum(sum(bb(ii).err));
end
% bb(end).out(2).val,out(2).val
err=err./prod(size(testimg));
twoerr=sort(err);

%% statistics
sprintf('%d, ',[aa.id,id]')
[~,idx]=min(err);
groundtruth=[idx,find(floor(([aa.id]-1)./10)==floor((id-1)/10))];
sprintf('%.2f, ',err')
% sparity concentration index
ii=idx;
out=aa(ii).out;
for jj=1:length(out),SCIval(jj)=norm(bb(ii).out(jj).val,1);end
% SCIval=(Nper*max(SCIval)/norm(out(jj).val,1)-1)/(Nper-1);
SCIval=(length(out)*max(SCIval)/norm(out(jj).val,1)-1)/(length(out)-1);
sprintf('%.2f, ',[SCIval,min(err),twoerr(2)/twoerr(1)])

%% display test recognition result
subplot(221),imshow(testimg,[]),title(sprintf('%d, ',groundtruth));
[X,Y]=meshgrid(1:4,1:4);
pos=[reshape(X,[],1) reshape(Y,[],1)];
orig=[];
for ii=1:4,tmp=[];
for jj=1:4,tmp=[tmp,aa((ii-1)*4+jj).im];end;orig=[orig;tmp];end
subplot(222),imshow(orig,[]);
orig=[];
for ii=1:4,tmp=[];
for jj=1:4,tmp=[tmp,bb((ii-1)*4+jj).im];end;orig=[orig;tmp];end
subplot(223),imshow(orig,[]);
orig=[];
for ii=1:4,tmp=[];
for jj=1:4,tmp=[tmp,bb((ii-1)*4+jj).err];end;orig=[orig;tmp];end
subplot(224),bar(err),% ;imshow(orig,[]),
title(sprintf('%.2f, ',[min(err),twoerr(2)/twoerr(1)]));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
drawnow;pause(2);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

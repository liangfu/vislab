function hogtrain_kfd
% addpath('traning');

system('../bin/hogtrain');

fp=fopen('../data/open.bin','r');
dat=fread(fp,10000000,'float');X=reshape(dat,378,[]);
Xtrain=X;N0=size(X,2)
fclose(fp);

fp=fopen('../data/close.bin','r');
dat=fread(fp,10000000,'float');X=reshape(dat,378,[]);
Xtrain=[Xtrain,X];N1=size(X,2)
fclose(fp);

Xtrain=Xtrain';

Y=zeros(size(Xtrain,1),1);
Y(N0+1:end)=1;

PRIOR=GPClearn_old(Xtrain,Y);

Wstr=sprintf('%ff,',PRIOR.alpha');Wstr=Wstr(1:end-1);
Pstr=sprintf('%ff,',PRIOR.X');Pstr=Pstr(1:end-1);

fp=fopen('../include/cvkfd4hog_template.h','r');
template=fread(fp,10000000,'char=>char');
fclose(fp);
out=sprintf(template,size(Xtrain,2),size(Xtrain,1),Wstr,Pstr);

fp=fopen('../include/cvkfd4hog.h','w');
fprintf(fp,out);
fclose(fp);


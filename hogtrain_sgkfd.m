function hogtrain_sgkfd
% addpath('traning');

system('../bin/hogtrain')

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

PRIOR=SGKFDlearn(Xtrain,Y); 

assert(size(PRIOR.alpha,2)==1);
Wstr=sprintf('%ff,',PRIOR.alpha');Wstr=Wstr(1:end-1); 
Pstr=sprintf('%ff,',PRIOR.X');Pstr=Pstr(1:end-1);

fp=fopen('../include/cvsgkfd4hog_template.h','r');
template=fread(fp,10000000,'char=>char');
fclose(fp);
out=sprintf(template,size(PRIOR.X,2),size(PRIOR.X,1),Wstr,Pstr,PRIOR.b);

fp=fopen('../include/cvsgkfd4hog.h','w');
fprintf(fp,out);
fclose(fp);


function hogtrain_svm

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

Y=-ones(size(Xtrain,1),1); 
Y(N0+1:end)=1;

fp=fopen('../dependency/libsvm/hogtraindata','w');
for ii=1:size(Xtrain,1)
datastr=[];
if Y(ii)>0,datastr=[datastr sprintf('+1 ')];
else,datastr=[datastr sprintf('-1 ')];end
for jj=1:size(Xtrain,2)
datastr=[datastr sprintf('%d:%f ',jj,Xtrain(ii,jj))];
end
datastr=[datastr sprintf('\n')];
fprintf(fp,'%s',datastr);
% disp(sprintf('%d,',ii));
end
fclose(fp);


% system('../dependency/libsvm/svm-train hogtraindata');
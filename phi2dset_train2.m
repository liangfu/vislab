function phi2dset_train

clear all;

nr=16;nc=16;
wdir='../data/palmdata/collection/all/';
files = dir([wdir '*.png']);
dat=[];
for ii=1:length(files)
  bw=imread([wdir files(ii).name])>0;
  im_rst = uint8(imnormalize(bw,size(bw,1)*size(bw,2)*0.25));
  phi = double(initsdf(im_rst));
  phi_dct=dct2(phi);
  phi_dct=phi_dct(1:nr,1:nc);
  dat(:,ii)=reshape(idct2(phi_dct)',[],1);
end
N=size(dat,2);


% X=dat';
for ii=1:N
tmp=reshape(dat(:,ii),nc,nr)';
tmp=tmp';
X(ii,:)=reshape(tmp,1,[]);
end
% imshow(reshape(mean(X),16,16)',[]);
% X(1,1:5)

% nframes=8;
% for ii=1:N-(nframes-1),
% stip(ii,:)=reshape(X(ii:ii+(nframes-1),:),1,[]);
% end

Y=zeros([N,1]);
% Y(329:522)=1;Y(562:589)=1;Y(656:579)=1;Y(796:843)=1;Y(851:906)=1;
% Y(1:5)=-1;Y(97:104)=-1;Y(491:501)=-1;Y(565:571)=-1;Y(845:850)=-1;
Y(1179:end)=1;

% trnidx=find(Y>=0);
% X=X(trnidx,:);Y=Y(trnidx);

% [W,PRIOR]=kfd(X(100:end-20,:),Y(100:end-20,:));
[W,PRIOR]=lda(X,Y,5);

hdrstr=['#ifndef __CV_CLASSIFIER_FOR_LEVELSETS_H__\n'...
        '#define __CV_CLASSIFIER_FOR_LEVELSETS_H__\n\n'...
        '#include "cvlda.h"\n\n'...
        'class CV_EXPORTS CvClassifier4Levelsets : public CvLDA'...
        '{\n'...
        ' public:\n'...
        '    CvClassifier4Levelsets():\n'...
        '  CvLDA()\n'...
        '  {\n'...
        '    static float W_data[]={\n'];
assert(size(W,1)==5&&size(W,2)==256);
Wstr=sprintf('%ff,',W');
hdrstr=[hdrstr Wstr(1:end-1) '\n'];
hdrstr=[hdrstr '    };\n\n'...
        '    static float PRIOR_data[]={\n'];
assert(size(PRIOR,1)==2&&size(PRIOR,2)==5);
PRIORstr=sprintf('%ff,',PRIOR');
hdrstr=[hdrstr PRIORstr(1:end-1) '\n};\n\n'...
        '    int Wnr=5;int Wnc=256;\n'...
        '    int PRIORnr=2;int PRIORnc=Wnr;\n\n'...
        '    assert(!W);assert(!PRIOR);\n'...
        '    W=cvCreateMat(5,256,CV_32F);\n'...
        '    PRIOR=cvCreateMat(2,5,CV_32F);\n'...
        '    memcpy(W->data.fl,W_data,sizeof(float)*Wnr*Wnc);\n'...
        '    memcpy(PRIOR->data.fl,PRIOR_data,'...
        '    sizeof(float)*PRIORnr*PRIORnc);\n'...
        '  }\n'...
        '};\n\n'...
        '#endif // __CV_CLASSIFIER_FOR_LEVELSETS_H__\n\n'];
fp=fopen('../include/cvclassifier4ls.h','w');
fprintf(fp,hdrstr);
fclose(fp);

opt = kpca(X,int32(size(X,2)),0); proj=opt.proj;
% [~,~,~,proj]=pca(X);
% W=kfd(real(proj),Y);
% W=kfd(X,Y);
% W=lda(X,Y);sprintf('%.5f,',W)
W=lda(real(proj),Y);
% W=lda_new(real(opt.proj),Y);

startp=1+49*4;
clf
for ii=startp:min(N,startp+48)
subplot(7,7,ii-startp+1),
imshow(reshape(X(ii,:),nc,nr)'>0),
title(num2str(ii))%,drawnow
end

% size(dat)
% mu(:,:)=mean(dat,3); 
% mu=mu';
% imagesc(mu);
% fclose(fp);


function hogtrain_lda

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

[W,PRIOR]=lda(Xtrain,Y,2);

Wstr=sprintf('%ff,',W');Wstr=Wstr(1:end-1);
Pstr=sprintf('%ff,',PRIOR');Pstr=Pstr(1:end-1);

out=[];
out=[out '#ifndef __CV_LDA4HOG_H__\n'];
out=[out '#define __CV_LDA4HOG_H__\n\n'];
out=[out '#include "cvlda.h"\n\n'];
out=[out 'class CV_EXPORTS CvLDA4HOG : public CvLDA\n'];
out=[out '{\n'];
out=[out 'public:\n'];
out=[out '  CvLDA4HOG():CvLDA()\n'];
out=[out '  {\n'];
out=[out '    static float W_data[]={\n'];
out=[out '      ' Wstr '\n'];
out=[out '    };\n'];
out=[out '    static float PRIOR_data[]={\n'];
out=[out '      ' Pstr '\n'];
out=[out '    };\n\n'];
out=[out '    int Wnr=2;int Wnc=378;\n'];
out=[out '    int PRIORnr=2;int PRIORnc=Wnr;\n\n'];
out=[out '    assert(!W);assert(!PRIOR);\n'];
out=[out '    W     = cvCreateMat(Wnr,Wnc,CV_32F);\n'];
out=[out '    PRIOR = cvCreateMat(PRIORnr,PRIORnc,CV_32F);\n\n'];
out=[out '    memcpy(W->data.fl,W_data,sizeof(float)*Wnr*Wnc);\n'];
out=[out '    memcpy(PRIOR->data.fl,PRIOR_data,sizeof(float)*PRIORnr*PRIORnc);\n'];
out=[out '  }\n'];
out=[out '};\n\n'];
out=[out '#endif // __CV_LDA4HOG_H__\n\n'];

fp=fopen('../include/cvlda4hog.h','w');
fprintf(fp,out);
fclose(fp);


function datwrite(X,fn)
% DATWRITE function writes matrix to .DAT file in the following
% format:
%    float 2 7291 256
%    <binary>
% Example:
%    X=randn([100,100]);X(20:80,20:80)=1;
%    datwrite(X,'../data/test.dat');
fp=fopen(fn,'w');
fprintf(fp,'%s %d ',...
        class(X),size(size(X),2));
fprintf(fp,'%d ',size(X));
fprintf(fp,'\n');
fwrite(fp,X',class(X));
fclose(fp);
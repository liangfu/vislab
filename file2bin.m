function file2bin

fp=fopen('test.txt','r');
data=[];N=-1;
while length(data)~=N,N=length(data);data=[data,fread(fp,1,'uchar')];end 
fclose(fp);
sprintf('0x%x,',data)
sprintf('%d,',data)

function partmatching

clear all; clf;

wdir='../dataset/yalefaces/';
ddir='';
gtruthfn='annotation.txt';

fp = fopen([wdir gtruthfn],'r');
gtruth=[];
I=1;N=0;
while 1
  gtruth(I).name = fscanf(fp, '%s',1); 
  if length(gtruth(I).name)==0 || gtruth(I).name(1)=='-', N=I-1;break; end
  gtruth(I).data = fscanf(fp, '%f',10);
  gtruth(I).data = reshape(gtruth(I).data,2,[])';
  I=I+1;
end
fclose(fp);

sz=32;
for I=1:N
data=gtruth(I).data;
scale=sqrt((data(1,1)-data(2,1)).^2+(data(1,2)-data(2,2)).^2)/(sz*1.5);
img = imread([wdir ddir gtruth(I).name]);
nr=size(img,1);nc=size(img,2);
clf;imshow(img,[]);
% hold on,plot([100,100],[100,200],'r-'),hold off;
npts=200;
hold on,plot(randi([1,nc],[npts,1]),randi([1,nr],[npts,1]),'r+'),hold off;
break;
end


function usps

clear all;

% fp=fopen('../data/usps','r');
% N=-1;X=[];Y=[];
% while 1
% N=size(X,1);
% t1=fscanf(fp,'%d ',1);
% t2=fscanf(fp,'%d:%f ',[2,256]);
% if length(t1)~=0||size(t2,2)~=0,X(N+1,:)=t2(2,:);Y(N+1,1)=t1;else,break;end
% end
% fclose(fp);
% size(X),size(Y)

datwrite([single(Y),single(X)],'usps.dat');
tmp=datread('usps.dat');
X=tmp(:,2:end);
Y=tmp(:,1);
size(X),size(Y)


C=unique(Y);
K=length(C);

% display
clf
for ii=1:K
subplot(2,5,ii),imshow(reshape(mean(X(Y==ii,:)),16,[])',[]);
end

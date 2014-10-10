function prob=hogdetect(I,W,PRIOR,csize,bound)
% HOGDETECT
% Usage :
%    prob=hogdetect(I,W,PRIOR,csize,bound)

if size(I,3)==3,I=rgb2gray(I);end

N=8;
%I=imresize(I,0.8);
[H,G]=hog(I,csize,bound);

nr=floor((size(I,1)-2)/csize);
nc=floor((size(I,2)-2)/csize);
assert(nr*nc*9==size(H,2))

cont=1;
for yy=1:nr-N
for xx=1:nc-N
for m=1:N
  Htest(cont,(m-1)*N*9+1:(m-1)*N*9+N*9)=...
      H(1,nc*9*(yy+m-1)+xx*9+1:nc*9*(yy+m-1)+xx*9+N*9);
end
cont=cont+1;
end
end

L=Htest*W';
%L=X(1:10,:)*W';
P1=exp(-.5*sum((L-repmat(PRIOR(1,:),size(L,1),1)).^2,2));
P2=exp(-.5*sum((L-repmat(PRIOR(2,:),size(L,1),1)).^2,2));
prob=zeros(size(L,1),2);
prob(:,1)=P1./(P1+P2);
prob(:,2)=P2./(P1+P2);
idx=find(prob(:,1)>0.66);

N=size(idx,1),
loc=[floor(idx/nc) mod(idx,nc)]
candidate=prob(idx)

clf,imagesc(G);


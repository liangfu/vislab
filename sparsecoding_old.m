function sparsecoding_old

if 1
  % load image
  X=double(rgb2gray(imread('../data/lena.jpg')))./255.0;
  X=imresize(X,1/2);
  X=double(X);

  Dsz=16; % dictionary size

  % build dictionary
  D=zeros([Dsz*Dsz,int32(size(X,1)/Dsz)*int32(size(X,1)/Dsz)]);
  nr=int32(size(X,1)/Dsz);nc=int32(size(X,2)/Dsz);
  for r=1:nr,
  for c=1:nc,
    D(:,(r-1)*nc+c)=reshape(X((r-1)*Dsz+1:(r-1)*Dsz+Dsz,...
                              (c-1)*Dsz+1:(c-1)*Dsz+Dsz),[Dsz*Dsz,1]);
  end
  end
  D_orig=D;
  D=D(:,randperm(size(D,2)));
end

% 
D=D_orig;
% [D,Di,beta,c,Q]=OMP(double(D(:,139)'),double(D),1e-8,10,121);
[D,Di,k]=OMP(double(D(:,139)'),double(D),1e-8,10,121);
Di(1:10),Di(k)

% show dictionary
if 1
  Dsz=int32(sqrt(size(D,1)));
  nr=int32(sqrt(size(D,2)));
  nc=int32(size(D,2))/nr;
  I = zeros(nr*Dsz,nc*Dsz);
  for r=1:nr,
  for c=1:nc,
    I((r-1)*Dsz+1:(r-1)*Dsz+Dsz,(c-1)*Dsz+1:(c-1)*Dsz+Dsz)=...
        reshape(D(:,(r-1)*nc+c),[Dsz,Dsz]);
  end
  end
  imshow(I,[]);
end



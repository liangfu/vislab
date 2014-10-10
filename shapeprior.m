function shapeprior

% pwdir = '../data/five/';
% files = dir([pwdir '*.pgm']);
pwdir = '../data/palmdata/collection/all/';
files = dir([pwdir '*.png']);

N = size(files,1);
M = [16 16];
objsz = 640;

im = imread([pwdir files(1).name]);
sz = size(im);
nr = sz(1);nc = sz(2);

im_sum = uint8(zeros(nr,nc));

for i = 1:N                             % for each training image

im = imread([pwdir files(i).name]);

%% normalize shape into common space
im_rst = uint8(imnormalize(im,objsz));
im_sum=im_sum+im_rst;

%% write normalized image data
%imwrite(im_rst.*255, [pwdir files(i).name '-norm.pgm']);

phi = double(initsdf(im_rst));
% rawwrite(['data/' files(i).name '-sdf.raw'], phi);

%% collect raw shape
shapeset(i,:,:)=phi(:,:);

end % end of scaning training dataset

%% compute mean shape 
tmp = mean(shapeset,1);
meanphi(:,:)=tmp(1,:,:);
%size(meanphi)
%imshow(meanphi,[]);

%% convert to DCT feature space
shapeVset=[];
for i=1:N
  tmp = [];tmp(:,:)=shapeset(i,:,:);
  shapeV = dct2(tmp'); %-meanphi
  shapeVtmp = shapeV(1:M(1),1:M(2));
  shapeVset(i,:) = reshape(shapeVtmp,1,[]);
end

%% compute mean shape
meanshapeV=mean(shapeVset,1);
featvec = reshape(meanshapeV,M(1),[]);
meanshapephi = ...
    idct2([featvec, zeros([M(1) nc-M(2)]); zeros([nr-M(1),nc])])';

%% perform PCA upon training set
[pc latent mu covar] = pca(shapeVset);
pc = pc(:,1:min(size(pc,2),N)); % dimension reduction
1-latent(1:20)./sum(latent)

%size(find(meanshapephi>0),1)
%gravitycenter(meanshapephi>0)

if 0 % iter==1, % save shape
  rawwrite(['../data/' 'shapeprior-meanshape.raw'], meanshapephi);
  rawwrite(['../data/' 'shapeprior-mean.raw'], meanshapeV); 
  % rawwrite(['../data/' 'shapeprior-covar.raw'], covar);
  rawwrite(['../data/' 'shapeprior-pc.raw'], pc(:,1:6)); % szpc=size(pc)
  rawwrite(['../data/' 'shapeprior-latent.raw'], latent); 
  % szlatent=size(latent)
else
  thedatastr=['/**\n * @file   cvshapeprior_data.cpp\n * @author ' ...
              'Liangfu Chen <liangfu.chen@cn.fix8.com>\n * @date   ' ...
              'Mon May 20 14:25:23 2013\n * \n * @brief  \n * \n * ' ...
              '\n */\n\n#include "cvshapeprior.h"\n\nvoid ' ...
              'CvShapePriorData::initialize()\n{\n  '...
              'assert(!meanshape);\n  ' ...
              'assert(!mean);\n  assert(!pc);\n  '...
              '// assert(!latent);\n\n  ' ...
              'float meanshape50x50_data[2500]={\n'];
  thedatastr=[thedatastr sprintf('%ff,',meanshapephi')];
  thedatastr=[thedatastr(1:end-1) '};\n\n  float mean256x1_data[256]={\n'];
  thedatastr=[thedatastr sprintf('%ff,',meanshapeV')];
  thedatastr=[thedatastr(1:end-1) '};\n\n  float pc12x256_data[3072]={\n'];
  thedatastr=[thedatastr sprintf('%ff,',pc(:,1:12)')]; 
  thedatastr=[thedatastr(1:end-1) '};\n\n  '...
'meanshape=cvCreateMat(50,50,CV_32F);\n  '...
'memcpy(meanshape->data.fl,meanshape50x50_data,sizeof(float)*50*50);\n  '...
'mean=cvCreateMat(1,256,CV_32F);\n  '...
'memcpy(mean->data.fl,mean256x1_data,sizeof(float)*256);\n\n  '...
'pc=cvCreateMat(256,12,CV_32F);\n  '...
'memcpy(pc->data.fl,pc12x256_data,sizeof(float)*12*256);\n}\n'];
  fp=fopen('../src/cvshapeprior_data.cpp','w');
assert(sum(size(meanshapephi)==[50,50])==2);
assert(sum(size(meanshapeV)==[1,256])==2);
assert(sum(size(pc(:,1:12))==[256,12])==2);
  fprintf(fp,thedatastr);
  fclose(fp);
  return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% shape similarities

% for i=1:N
% for j=1:N
%   similarVset(i,j)=1/(1+imdist(shapeVset(i,:), shapeVset(j,:)));
% end
% end
% similarVset=similarVset-diag(diag(similarVset));
% similarVset = similarVset./repmat(sum(similarVset,1),[N 1]);
% [junk idx] = max(similarVset,[],1); % idx

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% capture variation of a sampled shape

pc = pc(:,1:32); % dimension reduction

projset=[];
%for sampleid = 48:48
sampleid=48;
  
%% reconstruction
shapeV = reshape(shapeVset(sampleid,:),M(1),[]);
newshapeV = [ [shapeV, zeros([M(1) nc-M(2)])] ; zeros([nr-M(1),nc]) ];
newphi = idct2(newshapeV)';
%kappa = rbf(shapeVset(sampleid,:),meanshapeV,0);

%% variance ratio
ratio = cumsum(latent)/sum(latent); 
%ridx=find(ratio<1.0); ratio=ratio(ridx)
proj = (shapeVset(sampleid,:)-meanshapeV)*pc; % szproj=size(proj)
projset(sampleid,:)=proj;

%% reconstruction with lowest variance
X_var = proj*pc'+meanshapeV;
featvec = reshape(X_var, M(1), []);
varshapephi = idct2([[featvec,zeros([M(1), nc-M(2)])] ; ...
                    zeros([nr-M(1),nc])])';

%% calculate shape similarity
rms_error = imdist(varshapephi,newphi);

%% display
clf;
if 1
  subplot(151),imshow(del2(double(newphi<0))>0,[]),title('reconst. shape');
  subplot(152),imshow(double(newphi),[]),title('reconstructed \phi');
  subplot(153),imshow(real(varshapephi),[]),
  title(['pca #' num2str(size(pc,2)) ' reproj']);
  subplot(154),imshow(sqrt((varshapephi-newphi).^2),[]), ...
      title(['pca error: ' num2str(rms_error)]);
  subplot(155),imshow(meanshapephi,[]),title('mean \phi');
  % colormap('jet'); drawnow; pause(1);
  colormap('gray'); % drawnow; pause(1);
end

%end % end of enumerating samples

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% perform LDA upon low dimensional feature vectors 

% N1=7;N2=16;
% X=projset;
% W=lda(X,[zeros(N1,1); ones(N2,1)]);

% L = [ones(N1+N2,1) projset] * W';
% P = exp(L) ./ repmat(sum(exp(L),2),[1 2]);

% % ratio0=sum(P(1:N1,1)>0.5)/double(N1),
% % ratio1=(N2-sum(P(N1+1:N1+N2,1)>0.5))/double(N2)

% clf;
% hold on;
% h=plot(X(1:N1,1),X(1:N1,2),'rs');
% set(h,'MarkerFaceColor','r');
% h=plot(X(N1+1:N1+N2,1),X(1+N1:N1+N2,2),'bd');
% set(h,'MarkerFaceColor','b');
% hold off;


function sss

pwdir = '../data/five/';

files = dir([pwdir '*.pgm']);
N = size(files,1);
M = [16 32];

im = imread([pwdir files(1).name]);
sz = size(im);
nr = sz(1);nc = sz(2);

im_sum = uint8(zeros(nr,nc));

for i = 1:N                             % for each training image

im = imread([pwdir files(i).name]);

%% normalize shape into common space
im_rst = uint8(imnormalize(im,450.));
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
size(meanphi)
imshow(meanphi,[]);

%% convert to DCT feature space
for i=1:N
  tmp = [];tmp(:,:)=shapeset(i,:,:);
  shapeV = dct2(tmp-meanphi); %-meanphi
  shapeVtmp = shapeV(1:M(1),1:M(2));
  shapeVset(i,:) = reshape(shapeVtmp,1,[]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% shape similarities

for i=1:N
for j=1:N
  similarVset(i,j)=1/(1+imdist(shapeVset(i,:), shapeVset(j,:)));
end
end

similarVset=similarVset-diag(diag(similarVset));
similarVset = similarVset./repmat(sum(similarVset,1),[N 1]);
[junk idx] = max(similarVset,[],1); idx
subplot(155),imagesc(similarVset),colormap('gray'); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% capture variation of a sampled shape

sampleid = 4;

% reconstruction
shapeV = reshape(shapeVset(sampleid,:),M(1),[]);
size(shapeVset)
newshapeV = [ [shapeV, zeros([M(1) nc-M(2)])] ; zeros([nr-M(1),nc]) ];
newphi = idct2(newshapeV)+meanphi;

% compute mean shape
meanshapeV=mean(shapeVset,1);
featvec = reshape(meanshapeV,M(1),[]);
meanshapephi = ...
    idct2([featvec, zeros([M(1) nc-M(2)]); zeros([nr-M(1),nc])])+meanphi;
%kappa = rbf(shapeVset(sampleid,:),meanshapeV,0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for iter=1:10

% project to original data
[pc latent mu covar] = pca(shapeVset);
%size(find(meanshapephi>0),1)
%gravitycenter(meanshapephi>0)
if iter==1, % save shape
  rawwrite(['../data/' 'shapeprior-mean.txt'], meanshapephi);
  rawwrite(['../data/' 'shapeprior-covar.txt'], covar);
  rawwrite(['../data/' 'shapeprior-pc.txt'], pc);
  rawwrite(['../data/' 'shapeprior-latent.txt'], latent); 
end

pc = pc(:,1:iter); % dimension reduction
ratio = cumsum(latent)/sum(latent); % ratio(1:7)'
repmean = repmat(meanshapeV,[N 1]);
%proj = (shapeVset-repmean)*pc;
proj = (shapeVset(sampleid,:)-meanshapeV)*pc;

X_var = proj*pc'+meanshapeV;
%X_var = proj*pc'+repmean;
featvec = reshape(X_var, M(1), []);
varshapephi = idct2([[featvec,zeros([M(1), nc-M(2)])] ; ...
                    zeros([nr-M(1),nc])])+meanphi;

% calculate shape similarity
rms_error(iter) = imdist(varshapephi,newphi);
if rms_error(iter)<0.3, rms_error',break; end
% sum(sum( (abs(varshapephi)-abs(newphi)).^2 ))/(nc*nr)

% display
if 1
  subplot(151),imshow(del2(double(newphi<0))>0,[]),title('reconst. shape');
  subplot(152),imshow(double(newphi),[]),title('reconstructed \phi');
  subplot(153),imshow(real(varshapephi),[]),title('pca backproj');
  %  subplot(154),imshow(sqrt((varshapephi-newphi).^2),[]), ...
  %    title(['pca error: ' num2str(rms_error(iter))]);
  subplot(155),imshow(meanshapephi>0,[]),title('mean \phi');
  colormap('jet'); drawnow;
end

%if i==4,break;end
end




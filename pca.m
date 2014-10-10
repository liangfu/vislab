function [pc,latent,mu,proj] = pca(X)
% PCA - Perform PCA using covariance.
% data - MxN matrix of input data
% (M dimensions, N trials)
% signals - MxN matrix of projected data
% pc - each column is a PC
% latent - Mx1 matrix of variances
% Example : 
%   N=50;M=4;X=randn(N,M);
%   [pc,latent,mu,covar]=pca(X);
% or,
%   N=50;M=3;X=[rand(N,M)-0.5;rand(N,M)+0.5];[pc,latent,mu,covar]=pca(X);
[N M] = size(X);
% subtract off the mean for each dimension
mu = mean(X,1);
X_c = X-repmat(mu,[N 1]);
% calculate the covariance matrix
covar = cov(X); % covar = 1/N*X_c'*X_c;
% find the eigenvectors and eigenvalues
[pc, latent] = eig(covar');
% extract diagonal of matrix as vector
latent = diag(latent);
% sort the variances in decreasing order
[junk, idx] = sort(-1*latent);
latent = latent(idx);
pc = pc(:,idx);

if 1
% project the original data set
pc=pc(:,1:M);
proj = X_c*pc; % szproj = size(proj),proj(1,:),X(1,:);
% plot(proj(:,1),proj(:,2),'g.');
% back projection
X_var = proj*pc'+repmat(mu,[N 1]);
if 0,
  plot(proj(:,1),proj(:,2),'g.', ...
       X(:,1),X(:,2),'r.', ...
       X_var(:,1),X_var(:,2),'b.');
end
end
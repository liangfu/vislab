%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% coefs = get_spline_contours(xc)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function finds the B-spline coefficients for each of the
% object contours described by the contour parameters for each
% sample.  
%
% Input:
% xc: The contour parameters for each sample in the current frame. 
%
% Global variables required:
% k: the number of objects. 
% contour2coefs: The B-spline coefficients of the hand-labeled
% training contours.  contour2coefs is a 2 x (L + 3) x ncontours
% matrix s.t. contour2coefs(:,:,i) are the coefficients for contour
% number i.  
%
% Output:
% coefs: The B-spline coefficients of each object of each sample.  
% coefs(:,:,i,m) are the coefficients for object m of sample i.  
%
% This function extracts the affine transformation parameters from
% the contour parameter vector, then calls transform spline to
% apply this transformation to the training contour coefficients.  
%
% Functions called by this function:
% set_contour_indices
% transform_coefs
%
function coefs = get_spline_contours(xc)

set_contour_indices;
k = size(xc,1)/NPARAMS_PER_OBJECT;
global contour2coefs;

% number of parameters per object
nparams = size(xc,1)/k;

% number of samples
N = size(xc,2);

% Which element of xc is which
set_contour_indices;

% For each object
for m = 1:k,
  
  istart = (m-1)*nparams;
  
  % Find the transformations of the original contour
  mu = xc(istart+[CXLOC,CYLOC],:);
  a = xc(istart+CSMAJ,:);
  b = xc(istart+CSMIN,:);
  theta = (xc(istart+CTHETA,:) + pi * xc(istart+CISROT,:));
  isflip = xc(istart+CISFLIP,:);
  phi = xc(istart+CPHI,:);

  % Find the original contours
  contour = xc(istart+CCONTOUR,:);
  coefs0 = contour2coefs(:,:,contour);
  % Transform the original contours
  coefs(:,:,:,m) = transform_coefs(coefs0,mu,a,b,theta,phi,isflip);
  
end;

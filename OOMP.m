function [ D, Di, beta, c, Q ] = OOMP( f, D, tol, No, ind )
% OOMP  Optimized Orthogonal Matching Pursuit
%
% It creates an atomic decomposition of a signal using OOMP method [1]. You can choose a
% tolerance, the number of atoms to take in or an initial subspace to influence the OOMP
% algorithm. Non-selected atoms subtracted by their component in the chosen space are also
% available.
%  
% Usage:    [ Dnew, beta, Di ] = OOMP( f, D, tol );
%           [ Dnew, beta, Di, c, Q ] = OOMP( f, D, tol, No, ind ); 
%
% Inputs:
%   f       signal to be represented 
%   D       dictionary of atoms
%   tol     desired distance between f and its approximation the routine will stop if
%           norm(f'-Dsub*(f*beta)')*sqrt(delta)<tol where delta=1/L, L is number of points
%           in a sample
%   No      (optional) maximal number of atoms to choose, if the number of chosen atoms
%           equals to No, routine will stop (default No=size(D,2))
%   ind     (optional) indices determining  the initial subspace, 
%
% Outputs:
%   D       the dictionary D rearranged according to the selection process D(:,1:k)
%           contains the atoms chosen into the atomic decomposition 
%   beta    'k' biorthogonal functions corresponding to new D(:,1:k)
%   Di      indices of atoms in new D written w.r.t the original D
%   c       'k' coefficients of the atomic decomposition
%   Q       Q(:,1:k) contains orthonormal functions spanning new D(:,1:k), Q(:,k+1:N)
%           contains new D(:,k+1:N) subtracted by the projection onto the space generated
%           by  Q(:,1:k)  (resp. D(:,1:k))
%  
% References: 
%   [1] L. Rebollo-Neira and D. Lowe, "Optimized Orthogonal Matching Pursuit Approach",
%   IEEE Signal Processing Letters, Vol(9,4), 137-140, (2002). 
%   For the current implementation:
%   [2] M. Andrle and L. Rebollo-Neira, "A swapping-based refinement of orthogonal
%   matching pursuit strategies", Signal Processing, Vol 86, No 3, pp. 480-495, (2006). 
%
% See also OMP Swapping OOMPKSwaps OMPKSwaps OOMPFinalRefi BOOMP OBOMP

% See   http://www.ncrg.aston.ac.uk/Projects/HNLApprox/
%       http://www.nonlinear-approx.info/
name='OOMP';  %name of routine
% get inputs and setup parameters
[L,N]=size(D);
beta=[];
Di=1:N;     % index vector of full-dictionary 
Q=D; 
%delta=1/L;  %uncomment in the case of analytical norm
delta=1;   %uncomment in the case of discrete norm
if nargin<5 ind=[];end
if nargin<4 No=N;end
if nargin<3 tol=0.01;end; 

numind=length(ind);

%atoms having smaller norm than tol1 are supposed be zero ones
tol1=1e-7; %1e-5
%threshold for coefficients
tol2=1e-10;   %0.0001  %1e-5

tic;
fprintf('\n%s is running for tol=%g, tol1=%g and tol2=%g.\n',name,tol,tol1,tol2);
fprintf('tol  -> required precision, distance ||f-f_approx||\n')
fprintf('tol1 -> atoms having smaller norm are supposed to be zero ones.\n');
fprintf('tol2 -> algorithm stops if max(|<f,q>|/||q||)<= tol2.\n');
%===============================
% Main algorithm: at kth iteration
%===============================
H=min(No,N); %maximal number of function in sub-dictionary
for k=1:H    
  %finding of maximal coefficient
  c=zeros(1,N);cc=zeros(1,N);
  if k<=numind 
    [test,q]=ismember(ind(k),Di);
    if test~=1 error('Demanded index %d is out of dictionary',ind(k));end
    c(k)=norm(Q(:,q)); 
    if c(k)<tol1
      error('Demanded atom with index %d is dependent on the others.',ind(k));
    end;
  else 
    for p=k:N, c(p)=norm(Q(:,p));
      if c(p)>tol1 cc(p)=abs(f*Q(:,p))/c(p);end;
    end
    [max_c,q]=max(cc);
    %stopping criterion (coefficient)
    if max_c<tol2 
      k=k-1;
      fprintf('%s stopped, max(|<f,q>|/||q||)<= tol2=%g.\n',name,tol2);
      break;
    end
  end  
  if q~=k
    Q(:,[k q])=Q(:,[q k]); % swap k-th and q-th columns
    D(:,[k q])=D(:,[q k]);
    Di([k q])=Di([q k]);
  end	
  %re-orthogonalization of Q(:,k)  w.r.t Q(:,1),..., Q(:,k-1)         
  if k>1
    for p=1:k-1
      Q(:,k)=Q(:,k)-(Q(:,p)'*Q(:,k))*Q(:,p);
    end	   
  end
  nork=norm(Q(:,k)); 
  Q(:,k)=Q(:,k)/nork; %normalization
  % compute biorthogonal functions beta from 1 to k-1
  if k>1
    beta=beta-Q(:,k)*(D(:,k)'*beta)/nork; 
  end	
  beta(:,k)=Q(:,k)/nork;          % kth biorthogonal function
       
  %orthogonalization of Q(:,k+1:n) wrt Q(:,k)
  if k==N, break; end    
  for p=k+1:N
     % orthogonalization of Q(:,p) wrt Q(:,k)
     Q(:,p)=Q(:,p)-(Q(:,k)'*Q(:,p))*Q(:,k);
  end
  %stopping criterion (distance)
  if (tol~= 0) & (norm(f'-D(:,1:k)*(f*beta)')*sqrt(delta) < tol) break;end;
end

c=f*beta;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
normf=norm(f'-D(:,1:k)*c')*sqrt(delta);
fprintf('From %d atoms in this dictionary %s has chosen %d atoms.\n',N,name,k);
fprintf('\n%s lasted %g seconds\n',name,toc); 
fprintf('\nNorm ||f-f_approx|| is %g. \n',normf);
% error testing
% orthogonality
ErrorTest(Q(:,1:k));
% biorthogonality
ErrorTest(D(:,1:k),beta);

%Last modification Laura REBOLLO-NEIRA (2009) (former OOMPF)
%
%Copyright (C) 2006 Miroslav ANDRLE and Laura REBOLLO-NEIRA
%
%This program is free software; you can redistribute it and/or modify it under the terms 
%of the GNU General Public License as published by the Free Software Foundation; either 
%version 2 of the License, or (at your option) any later version.
%
%This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
%without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
%See the GNU General Public License for more details.
%
%You should have received a copy of the GNU General Public License along with this program;
%if not, write to the Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
%Boston, MA  02110-1301, USA.

% function [ D, Di, beta, c, Q ] = OMP( f, D, tol, No, ind )
function [ D, Di, k ] = OMP( f, D, tol, No, ind )
% OMP Orthogonal Matching Pursuit
%
% It creates an atomic decomposition of a signal using OMP criterion. You can choose a
% tolerance, the number of atoms to take in or an initial subspace to influence the OMP
% algorithm.  
% 
%  
% Usage:    [ Dnew, Di, beta, c ] = Omp( f, D, tol, No, ind );
%           [ Dnew, Di ] = Omp( f, D, tol );
%
% Inputs:
%   f       analyzing signal
%   D       dictionary of normalized atoms
%   tol     desired distance between f and its approximation the routine will stop if
%           norm(f'-Dsub*(f*beta)')*sqrt(delta)<tol where delta=1/L, L is number of points
%           in a sample
%   No      (optional) maximal number of atoms to choose, if the number of chosen atoms
%           equals to No, routine will stop (default No=size(D,2) 
%   ind     (optional) indices determining  the initial subspace, 
%
% Outputs:
%   D       the dictionary D rearranged according to the selection process D(:,1:k)
%           contains the atoms chosen into the atomic decomposition 
%   Di      indices of atoms in new D written w.r.t the original D
%   beta    'k' biorthogonal functions corresponding to new D(:,1:k)
%   c       'k' coefficients of the atomic decomposition
%  
% References: 
%   L. Rebollo-Neira and D. Lowe, "Optimized Orthogonal Matching Pursuit Approach", IEEE
%   Signal Processing Letters, Vol(9,4), 137-140, (2002). 
%
% See also OMPF.

% See   http://www.ncrg.aston.ac.uk/Projects/HNLApprox/
%       http://www.nonlinear-approx.info/

  
name='OMP';  %name of routine
% get inputs and setup parameters
[L,N]=size(D);
zmax=2;%number of reorthogonalizations
beta=[];
Re=f;
Di=1:N;     % index vector of full-dictionary 
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

% fprintf('\n%s is running for tol=%g, tol1=%g and tol2=%g.\n',...
%         name,tol,tol1,tol2);
% fprintf('tol  -> required precision, distance ||f-f_approx||\n')
% fprintf('tol1 -> '...
%         'atoms having smaller norm are supposed to be zero ones.\n');
% fprintf('tol2 -> algorithm stops if max(|<f,q>|/||q||)<= tol2.\n');

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
  else 
    for p=k:N, %it can be done without for
      cc(p)=abs(Re*D(:,p));
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
   D(:,[k q])=D(:,[q k]); % swap k-th and q-th columns
   Di([k q])=Di([q k]);
  end
  if k>1
   %Q(:,k) is the orthogonalization of D(:,k)  w.r.t Q(:,1),..., Q(:,k-1) 
    for p=1:k-1
     Q(:,k)=D(:,k)-(Q(:,p)'*D(:,k))*Q(:,p); %orthogonalization
    end 
   %re-orthogonalization of Q(:,k)  w.r.t Q(:,1),..., Q(:,k-1) 
   for zi=1:zmax
    for p=1:k-1
     Q(:,k)=Q(:,k)-(Q(:,p)'*Q(:,k))*Q(:,p); %re-orthogonalization
    end	   
   end
  end
  if k==1 Q(:,k)=D(:,k); end
  nork=norm(Q(:,k)); 
  Q(:,k)=Q(:,k)/nork; %normalization
  % compute biorthogonal functions beta from 1 to k-1
  if k>1
    beta=beta-Q(:,k)*(D(:,k)'*beta)/nork; 
  end	
  beta(:,k)=Q(:,k)/nork; % kth biorthogonal function
  %
%%%%computes de residum%%%%%   
  Re=Re-f*Q(:,k)*Q(:,k)';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%
  %stopping criterion (distance)
  if (tol~= 0) & (norm(f'-D(:,1:k)*(f*beta)')*sqrt(delta) < tol) break;end;
end

c=f*beta;
normf=norm(f'-D(:,1:k)*c')*sqrt(delta);

% fprintf('From %d atoms in this dictionary %s has chosen %d atoms.\n',...
%         N,name,k);
% fprintf('\n%s lasted %g seconds\n',name,toc); 
% fprintf('\nNorm ||f-f_approx|| is %g. \n',normf);

%% error testing
% orthogonality
% ErrorTest(Q(:,1:k));
% biorthogonality
% ErrorTest(D(:,1:k),beta);

%Last modification  Laura REBOLLO-NEIRA (2009)
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

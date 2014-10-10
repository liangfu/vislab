function [ D, Di, beta, c ] = BOOMP( f, D, beta, tol, No )
% BOOMP  Backward-Optimized Orthogonal Matching Pursuit
%
% Using the Least Square criterion at each step it eliminates one function from a given 
% basis to have best possible representation in the reduced space. It also modifies the
% corresponding biorthogonal functions.
%  
% Usage:    [ D, Di, beta, c ] = BOOMP( f, D, beta, tol, No );
%
% Inputs:  
%   f       signal to represent
%   D       set of chosen independent functions
%   beta    set of biorthogonal functions to D    
%   tol     tolerance, desired difference between signal and its approx, (optional,
%           default tol=1.0e-2)  
%   No      (optional) desired number of atoms in the decomposition if you want really No
%           atoms set tol=0, it speeds the process  
%  
% Outputs: 
%   D       new reduced set of independent functions
%   Di      indices of atoms in new D w.r.t. to original D
%   beta    biorthogonal functions to new D 
%   c       coefficients of the atomic decomposition  
%
% Note: this routine should be use only at the end of our selection process since it is
% not adapting (for the speed purposes) unselected dictionary functions. Thus any of our
% forward selection methods cannot be used after this.  
%  
% References:
%   M. Andrle, L. Rebollo-Neira, and E. Sagianos, "Backward-Optimized Orthogonal Matching 
%   Pursuit Approach", IEEE Signal Processing Letters, Vol (11,9), 705-708 (2004). 

% See   http://www.ncrg.aston.ac.uk/Projects/HNLApprox/
%       http://www.nonlinear-approx.info/
    
name='BOOMP';
tic;
[L,N]=size(D); 
if nargin==3; tol=1.0e-2;end; 
if nargin==4; No=0;end;
if No>N; error('%s: You want to reduce %d atoms into %d atoms!',name,N,No);end
%initilization
delta=1; %delta=1/L;
Di=1:N;
C=f*beta;c=C;   
nor1=norm(f'-D*C')*sqrt(delta);nor2=nor1;
if (nor1>tol) & (No==0) 
  fprintf('\n%s: Approximation is too coarse, I cannot satisfy the tolerance condition\n',name);
else
  for n=1:N-No
    s=size(beta,2);
    dd=zeros(1,s);%clear dd
    for i=1:s  
      dd(i)=beta(:,i)'*beta(:,i);
    end	
    [Cp,p]=min((C.*conj(C))./dd); % p-th atom to be removed
    temp=(beta'*beta(:,p))*C(p)/dd(p);
    C=C-temp';
    beta1=beta-beta(:,p)*(beta(:,p)'*beta)/dd(p);
    %stopping criterion
    if (tol~= 0) & (norm(f'-D*C')*sqrt(delta)>tol); break; end
    %if norm(f'-D*(f*beta1)')*sqrt(delta)>tol break; end
    beta=beta1;
    beta(:,p)=[];
    C(p)=[];c=C;
    D(:,p)=[];Di(p)=[];
  end
  nor2=norm(f'-D*c')*sqrt(delta);
end  
fprintf('\n%s reduces the space of %d functions into a space of %d functions\n',name,N,size(D,2));
fprintf('Previous norm ||f-f_approx|| was %g, new norm is %g\n',nor1,nor2);
if No~=0 
  fprintf('You have required %d atoms in decomposition\n',No);
else  
  fprintf('Required precision was tol=%g\n',tol);
end  
fprintf('%s ran for %f seconds\n',name,toc);

%Last Modification Laura REBOLLO-NEIRA (2009)
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

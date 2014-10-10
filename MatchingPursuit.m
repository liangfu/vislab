function MatchingPursuit
% MATCHINGPURSUIT - match pursuit demo
% Copyright (C) Liangfu Chen (2012-2013)

clear all;
D = [1/2 sqrt(3)/2; 1 0; -1/sqrt(2) -1/sqrt(2)]';
y = [1 1/2]';

%% estimation
err=y;
for iter=1:100
for ii=1:size(D,2),a(ii)=sum(err(:,iter).*D(:,ii));end
[maxval,maxloc]=max(a); 
out.val(iter)=maxval;out.loc(iter)=maxloc;
err2=err(:,iter)-maxval*D(:,maxloc);
err=[err,err2];
if(norm(err2)<.02),break;end
end

%% evaluation by reconstructing the signal
approx=zeros(2,1);
for ii=1:length(out.loc)
approx=approx+D(:,out.loc(ii))*out.val(ii);
end
[approx,y]

function phi0 = initsdf(bw)
% INITSDF - initialize level-set with signed distance function
if 1
  tmp = bwdist(bw);
  phi0 = bwdist(bw==0)-tmp;
else
  phi0=(bw>0)-(bw<1);
end
  
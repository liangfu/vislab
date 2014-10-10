function pdist = imdist(phi,phi0)
% DISTANCE - distance between levelsets
[nr nc]=size(phi);
pdist = sum(sum( abs(phi-phi0) ))/(nc*nr);

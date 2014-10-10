function H = Heaviside(phi, sigma)
  idx0 = find((phi<=sigma) & (phi>=-sigma));
  idx1 = find(phi<-sigma);
  idx2 = find(phi>sigma);
  H = zeros(size(phi,1),size(phi,2));
  H(idx0) = ((0.5/sigma)*phi(idx0) + 0.5/pi*sin(pi*phi(idx0)/sigma))+0.5; 
  H(idx1)=1e-5;                                  
  H(idx2)=1-1e-5;                                 


function gp

clear all;
xs = (-15:0.2:15)'; ns = size(xs,1); keps = 1e-9;
K = inline('exp(-0.5*(repmat(p'',size(q))-repmat(q,size(p''))).^2)');

% the unknown function to be estimated 
fs = cos(xs/pi)+randn(ns,1).*.2;

% generate testing points
x=(-14.3:.2:14.3)';
clf;
plot(xs,fs,'b.'); hold on;
for ii=1:length(x)
md=K(x(ii),xs)'*inv(K(xs,xs)+keps*eye(ns))*fs;
plot(x(ii),md(end),'r+'); 
end
hold off;


function H = gabor(opt, t, tau)
% GABOR - generate gabor kernel with given angle and scale
% Example:
%    opt='even';t=-15:0.1:15;tau=8;
%    h = gabor(opt, t, tau)

if nargin==0,opt='even';sz=4;tau=2;end  

omega=1/tau;
if strcmp(opt,'even')
  H=-cos(2*pi*t*omega).*exp(-(t.^2)/(tau^2));
elseif strcmp(opt,'odd')
  H=-sin(2*pi*t*omega).*exp(-(t.^2)/(tau^2));
end

if nargin==0,plot(t,H);end
if nargout==0,H=[];end
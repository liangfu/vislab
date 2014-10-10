function [im,opts,flog,plog,sclog] = superres_ml(W,Y,La,Lb,avim,opts)
%
% function [im,opts,flog,plog,sclog] = superres_ml(W,Y,La,Lb,avim,opts);
%
% W: W matrix (N*KM)
% Y: Y vector (KM*1)
% La: multiplicative photometric parameter.
% Lb: additive photometric parameter.
% avim: average image (biv*bih)
% opts: optional options vector to pass to scg.
%

GRADTEST = 0; % Set this to 1 in order to run local_gradcheck below.

if (nargin<6)
    fprintf('Using default options vector.\n');
    opts = zeros(1,18);
    opts(1) = 1;
    opts(2) = 1e-10;
    opts(3) = 1e-10;
    opts(14) = 40;
end
biv = size(avim,1);
bih = size(avim,2);

Y = Y-Lb;
KM = numel(La);
W = mex_amub(W,sparse(1:KM,1:KM,La));
% The line above would be "W = W*sparse(1:KM,1:KM,La);", but matlab doesn't
% handle memory very well in this case, so it generates a completely
% unnecessary "out of memory" error. I got mex_amub from the file exchange,
% though now "ssmult" seems to be preferred.

x = avim(:)';

if (GRADTEST)
    [D,G] = local_gradcheck(x,W,Y);
    im = D; opts = G; % overwrite usual outputs when testing grad.
    return;
end


if (opts(1)>=0), fprintf('hmem_ml iterating (%i iters):\n',opts(14)); end;
if (nargout <4)
    [X,opts,flog] = scg(@local_func,x,opts,@local_grad,W,Y);
else
    [X,opts,flog,plog,sclog] = scg(@local_func,x,opts,@local_grad,W,Y);
end
im = reshape(X,biv,bih);
if (opts(1)>=0), fprintf('Done!\n'); end;
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function e = local_func(x,W,Y)
e = Y'-x*W;
e = sum(e.^2);
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function g = local_grad(x,W,Y)
r = Y'-x*W;
g = W*r';
g = -2.*g';
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [D,G] = local_gradcheck(x,varargin)
% This function checks the gradient numerically, for a random subsection of
% the image pixels (since it takes rather too long to do all of them as in
% scg's automatic gradient test option, selected by setting opts(9)=1).

ep = 1e-4; % numerical grdient step size.
numsamp = 30; % approximate number of pixels I'd like to check.

x = x + reshape(gsamp(0,0.15^2,numel(x)),size(x)); % Add in a bit of randome noise.
G = local_grad(x,varargin{:});
L = numel(x);
if (L~=numel(G)), error('lengths of parameter vector and gradient vector are inconsistent.'); end;
step = zeros(1,L);
D = step;
t = numsamp/L;
for i = 1:L
    if (mod(i,10000)<1), fprintf('Element %i of %i.\n',i,L); end;
    if(rand<t)
        step(i) = ep;
        l1 = local_func(x+step,varargin{:});
        l2 = local_func(x-step,varargin{:});
        D(i) = (0.5/ep)*(l1-l2);
        step(i) = 0;
    else
        D(i) = NaN;
    end
end
fprintf('Done all elements!\n');
m = logical(~isnan(D));
G = G(m);
D = D(m);
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xi,yi]=gdsolver(fx,dfx,x0)
% GDSOLVER demostrates gradient descent algorithm on Wikipedia
    
if nargin<1
  fx{1}=@(x) 3*x(1)-cos(x(2)*x(3))-1.5;
  fx{2}=@(x) 4*x(1)*x(1)-625*x(2)*x(2)+2*x(2)-1;
  fx{3}=@(x) exp(-x(1)*x(2))+20*x(3)+(10*pi-3)/3.;
  dfx{1}=@(x) [3,sin(x(2)*x(3))*x(3),sin(x(2)*x(3))*x(2)];
  dfx{2}=@(x) [8*x(1),-1250*x(2)+2,0];
  dfx{3}=@(x) [-x(2)*exp(-x(1)*x(2)),-x(1)*exp(-x(1)*x(2)),20];
  x0=zeros([3,1]);
end

Gx=@(x,fx) (cellfun(@(c) c(x),fx))';      % `Gx(x0,fx)` - nonlinear system
Fx=@(x,Gx,fx) .5*Gx(x,fx).*Gx(x,fx);      % `Fx(x0,Gx,fx)` - objective function
J_G=@(x,dfx) cell2mat(cellfun(@(c) (c(x))',dfx,'uniformoutput',false)); % `J_G(x0,dfx)`
dFx=@(x,J_G,dfx,Gx,fx) (J_G(x,dfx))'*Gx(x,fx); % `dFx(x0,J_G,dfx,Gx,fx)` - derivative

% parameters
gamma=1e-3;
epsilon=0.1;
maxiter=1000;

% initialize
if ~exist('x0','var'),x0=zeros([length(fx),1]);end
xi=x0;
yi=sum(Fx(xi,Gx,fx)); % eval objective
iter=2;

% start iteration
while iter<maxiter && yi(end)>epsilon
xi=xi-gamma*dFx(xi,J_G,dfx,Gx,fx);
yi(iter)=sum(Fx(xi,Gx,fx)); % eval objective
iter=iter+1;
end
plot(yi);yi=yi(end);

end


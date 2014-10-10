function snake
% EXAMPLE     an example of traditional snake on U-shape image

[I,C]=testimg();
I=I*255;
k=fspecial('gaussian',[17,17],11);
I=imfilter(I,k);

% Compute its edge map, 
f = 1 - I/255; 
f0 = f; %gaussianBlur(f,1);
        % note: snake potential is the negative of edge map
[px,py] = gradient(f0);

% display the results
x=C(:,1);
y=C(:,2);
[x,y] = snakeinterp_v2(x,y);
clf,imshow(-f,[]);hold on,plot(x,y,'r-'),hold off;drawnow;

for i=1:100,
  [x,y] = snakedeform(x,y,0.05,0,1,4,px,py,5);
  [x,y] = snakeinterp_v2(x,y);

  clf,imshow(-f,[]);hold on,plot(x,y,'r-'),
  title(['Deformation in progress,  iter = ' num2str(i*5)]),hold off;drawnow;
  % break;
end

function [img,c]=testimg
% generate fixed sized test image
bw=zeros([256,256]);
nr=size(bw,1);
nc=size(bw,2);
noisy=int32(nc*nr*.3);
bw(int32(nr/4):int32(nr/4*3),int32(nc/4):int32(nc/4*3))=1;
x=randi([1,255],[noisy 2]);y=x(:,1)*256+x(:,2);
bw(y)=1;
x=randi([1,255],[noisy 2]);y=x(:,1)*256+x(:,2);
bw(y)=0;
img=double(bw);
% generate initial curve
angle=1e-5:.2:pi*2.;
tmp=randi([38,48],[1,max(size(angle))]);
c=int32([tmp.*cos(angle)+80;tmp.*sin(angle)*1.8+120]');
% display them, optional
imshow(img,[]);hold on,plot(c(:,1),c(:,2),'r-','LineWidth',2),hold off;
c=double(c);

function [x,y] = snakedeform(x,y,alpha,beta,gamma,kappa,fx,fy,ITER)
% SNAKEDEFORM  Deform snake in the given external force field
%     [x,y] = snakedeform(x,y,alpha,beta,gamma,kappa,fx,fy,ITER)
%
%     alpha:   elasticity parameter
%     beta:    rigidity parameter
%     gamma:   viscosity parameter
%     kappa:   external force weight
%     fx,fy:   external force field
% 
% Chenyang Xu and Jerry L. Prince, 4/1/95, 6/17/97
% Copyright (c) 1995-97 by Chenyang Xu and Jerry L. Prince
% Image Analysis and Communications Lab, Johns Hopkins University

% generates the parameters for snake

N = length(x);

alpha = alpha* ones(1,N); 
beta = beta*ones(1,N);

% produce the five diagnal vectors
alpham1 = [alpha(2:N) alpha(1)];
alphap1 = [alpha(N) alpha(1:N-1)];
betam1 = [beta(2:N) beta(1)];
betap1 = [beta(N) beta(1:N-1)];

a = betam1;
b = -alpha - 2*beta - 2*betam1;
c = alpha + alphap1 +betam1 + 4*beta + betap1;
d = -alphap1 - 2*beta - 2*betap1;
e = betap1;

% generate the parameters matrix
A = diag(a(1:N-2),-2) + diag(a(N-1:N),N-2);
A = A + diag(b(1:N-1),-1) + diag(b(N), N-1);
A = A + diag(c);
A = A + diag(d(1:N-1),1) + diag(d(N),-(N-1));
A = A + diag(e(1:N-2),2) + diag(e(N-1:N),-(N-2));

invAI = inv(A + gamma * diag(ones(1,N)));

for count = 1:ITER,
  vfx = interp2(fx,x,y);
  vfy = interp2(fy,x,y);
  
  % deform snake
  x = invAI * (gamma* x + kappa*vfx);
  y = invAI * (gamma* y + kappa*vfy);
end

function [xi,yi] = snakeinterp_v2(x,y)

x=x(:);y=y(:);
U=x+i*y;

% tranform
G=fft(U);

% normalize
center=G(1)/length(G);G(1)=0; % translate
scale=abs(G(2));G=G./scale; % scale

% modify
L_src=floor(length(G));
L_dst=floor(arclength([x,y])/1.5);
if L_dst>L_src
  Gi=[0;G(2:floor(L_src/2));zeros([L_dst-L_src,1]);G(end-floor(L_src/2-2):end)];
else
  Gi=[0;G(2:floor(L_dst/2));G(end-floor(L_dst/2-2):end)];
end

% reproduce
Gi=Gi.*scale*length(Gi)/length(G); % scale
Gi(1)=center*length(Gi); % translate

% transform
Ui=ifft(Gi);

xi=real(Ui);
yi=imag(Ui);

function L=arclength(c)
L=0;
for ii=1:size(c,1)-1
  L=L+sqrt(sum((c(ii,:)-c(ii+1,:)).^2));
end


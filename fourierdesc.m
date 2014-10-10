function fourierdesc

[I,C]=testimg();

x=double(C(:,1));
y=double(C(:,2));
U=x+i*y;

% tranform
G=fft(U);

% normalize
center=G(1)/length(G);G(1)=0; % translate
scale=abs(G(2));G=G./scale; % scale
theta=atan2(imag(G(2)),real(G(2))); % rotation
rx=real(G(:));ry=imag(G(:));
Gt=[cos(-theta) -sin(-theta);
    sin(-theta) cos(-theta)]*[rx';ry'];
G=Gt(1,:)+i*Gt(2,:);G=G(:);

% modify
L=floor(length(G)/2);
Gi=[0;G(2:L);zeros([200,1]);G(end-(L-2):end)];
%Gi=[0;G(2:L);G(end-L+2:end)];

% reproduce
angle=atan2(imag(Gi(2)),real(Gi(2))); % rotation
rx=real(Gi(:));ry=imag(Gi(:));
Gt=[cos(theta-angle) -sin(theta-angle);
    sin(theta-angle) cos(theta-angle)]*[rx';ry'];
Gi=Gt(1,:)+i*Gt(2,:);Gi=Gi(:);
Gi=Gi.*scale*length(Gi)/length(G); % scale
Gi(1)=center*length(Gi); % translate

% transform
Ui=ifft(Gi);

% display
clf;hold on;
plot(real(U ),imag(U ),'r-');
plot(real(Ui),imag(Ui),'r-');
hold off;

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


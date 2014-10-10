function chamfer

workdir = '../data/five/';
files = dir([workdir '*.pgm']);
N=length(files);

%for f = 4:N
f=4
im = imread([workdir files(f).name]);
%im=imnormalize(im,450);
im = double(im==0);
bw = del2(im)>0;
%subplot(3,3,f-3),
tmplt=zeros(121,121);
mask=zeros(121,121);
tmplt(35:35+50-1,35:35+50-1)=bw;
mask(35:35+50-2,35:35+50-1)=(im(2:end,:)==0);
imshow(tmplt);
imshow(mask);
rawwrite('../data/shapeprior-chamfer.raw',single(tmplt));
%rawwrite('../data/shapeprior-boundary.raw',single(tmplt));
%rawwrite('../data/shapeprior-region.raw',single(mask));

%break;

[X,Y]=meshgrid(-25:25,-25:25);
g = 5*exp(-(0.1*X.^2+0.03*Y.^2))-1.8*exp(-(0.005*X.^2+0.005*Y.^2))-0.5;
rawwrite('../data/shapeprior-gaussian.raw', g);
surf(X,Y,rawread('../data/shapeprior-gaussian.raw')),colormap('jet');

%break;

%end


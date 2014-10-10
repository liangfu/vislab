function randsample

f=dir('scene');
sz=[62,58];

N=5000;

out=[];outcc=1;
imgidx=randi([1,length(f)],[1 N]);
for ii=1:N
if f(imgidx(ii)).name(1)=='.',continue;end
im=rgb2gray(imread(['scene/' f(imgidx(ii)).name]));
nr=size(im,1);nc=size(im,2);
yy=randi(nr-sz(1)-2,1);
xx=randi(nc-sz(2)-2,1);
roi=im(yy:yy+sz(1)-1,xx:xx+sz(2)-1);
imwrite(roi,sprintf('neg/neg%04d.pgm',outcc));
out(outcc).roi=roi;outcc=outcc+1;
end


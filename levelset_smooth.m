function levelset_smooth
aa=5:-1:-5;bb=4:-1:-6;cc=6:-1:-4;
phi = [aa;aa;aa;aa;cc;aa;aa;aa;aa];
[dx dy]=gradient(phi);
f = div(dx./(sqrt(dx.^2+dy.^2)+1e-5), ...
        dy./(sqrt(dx.^2+dy.^2)+1e-5));
d = Dirac(phi,3);
lambda=5;
subplot(221),surf(phi);
subplot(222),imshow(phi+lambda*d.*f,[]);
subplot(223),imshow(lambda*d.*f,[]);
subplot(224),imshow((phi+lambda*d.*f)<0,[]);
colormap('jet');
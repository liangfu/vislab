function toydenoising

im = imread('../data/noisytoy.pgm');
bw=im>0;
m0=sum(bw,1);
m1=sum(bw,2);
subplot(311),imshow(bw,[]);
subplot(312),plot(m0);
subplot(313),plot(m1);


function dsnake

files = dir('data/*.pgm');
N = size(files,1);

im = imread(['data/' files(1).name]);
sz = size(im);
avg=zeros(sz);

for i = 1:N                             % for each training image

im = imread(['data/' files(i).name]);
sz = size(im);
bw = zeros(sz);

phi= initsdf(im); % rawwrite(['data/' files(i).name '-sdf.raw'], phi);
[yy xx] = find(phi<0);
for i = 1:size(yy)
  bw(yy(i),xx(i))=1;
end

ty=sz(1)/2-mean(yy);
tx=sz(2)/2-mean(xx);

for i = 1:sz(1)
for j = 1:sz(2)
  if ~((int32(i+ty)<1) || (int32(j+tx)<1) || ...
       (int32(i+ty)>sz(1)) || (int32(j+tx)>sz(2)))
    avg(i,j) = bw(int32(i+ty),int32(j+tx))+avg(i,j);
  end
end
end

end

subplot(1,2,1),imshow(avg,[]); 
subplot(1,2,2),imshow(avg,[]); 


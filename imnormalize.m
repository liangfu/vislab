function bw0 = imnormalize(bw, stdarea)
% IMNORMALIZE - normalize image to assgined area and gravity center
% bw0 = imnormalize(bw, stdarea)
% bw is the original input binary pattern
% stdarea is the desired size of the centered object

stdcy = size(bw,1)/2.; 
stdcx = size(bw,2)/2.;

% translate
idx0 = find(bw>0);
bw(idx0)=1;

% scaling
area = sum(sum(bw));
scale=sqrt(stdarea)/sqrt(area);
tmp=imresize(bw,scale,'cubic');

% adjust position
sz=size(bw);
newsz=size(tmp);
bw0=zeros(sz);
if scale>=1.
  szdiff=round((newsz-sz)./2.);
  bw0=tmp(szdiff(1)+1:szdiff(1)+sz(1), ...
          szdiff(2)+1:szdiff(2)+sz(2));
else
  szdiff=round((sz-newsz)./2.);
  bw0(szdiff(1)+1:szdiff(1)+newsz(1), ...
      szdiff(2)+1:szdiff(2)+newsz(2)) = tmp;
end

assert(size(bw,1)==size(bw0,1));
assert(size(bw,2)==size(bw0,2));
[center.y center.x]=imcentroid(bw0);
tform=maketform('affine', ...
                [1 0; 0 1; stdcx-center.x stdcy-center.y]);
bw0=imtransform(bw0,tform,'xdata',[1 sz(1)],'ydata',[1 sz(2)]);

% display
if 0
  [center.y center.x]=imcentroid(bw0);
  out.cy=center.y;out.cx=center.x;
  out.sumbw0=sum(sum(bw0));
  out
end



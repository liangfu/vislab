function bw=bwcenter(bw,sz,scale)
% BECENTER function centers a binary image pattern 
%    bw = bwcenter(bw,sz,scale)
% 
% Example :
%    bw = bwcenter(bw,sz,scale)

if nargin<2,if nargin<3,scale=0.25;end,end
bwsize = size(bw);
im=imnormalize([bw,zeros([bwsize(1),50-bwsize(2)]);...
                zeros([50-bwsize(1),50])],50*50*0.25);
bw=im>(max(max(im))/2.0);
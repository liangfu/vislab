function subimg=getsubrect(img,center,sz)
% GETSUBRECT sample ROI with given center and size
%
% Description:
%    SUBIMG=GETSUBRECT(IMG,CENTER,SZ) takes single channel
%    grayscale image IMG with center and size, generates defined
%    region of interest as output matrix.
% 
% Example:
%    img = imread('test.pgm');
%    center=[50,50];
%    sz=[10,10];
%    subimg=getsubrect(img,center,sz)
% 
% Copyright (c) Liangfu Chen (2012-2013)

subimg=img(int32(center(2)-sz(2)*.5):int32(center(2)+sz(2)*.5),...
           int32(center(1)-sz(1)*.5):int32(center(1)+sz(1)*.5));
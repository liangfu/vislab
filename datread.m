function X=datread(fn)
% DATREAD function reads matrix from .DAT file in the following
% format:
%    float 2 7291 256
%    <binary>
% Example:
%    X=datread('../data/test.dat');
fp=fopen(fn,'r');
type=fscanf(fp,'%s ',1);
dim =fscanf(fp,'%d ',1)
sz  =fscanf(fp,'%d ',dim)
X=fread(fp,prod(sz),type);
X=reshape(X,[],sz(1))';
fclose(fp);

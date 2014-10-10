function rawwrite(fn,img)
% RAWWRITE write a single matrix in the follow format:
%   nr nc
%   dat1 dat2 ...
% Example:
%   rawwrite("file.txt", single(img));
fid = fopen(fn,'w');
sz=size(img);
fprintf(fid,'%d %d\n',sz(1),sz(2));
fprintf(fid,'%f ',img');
fclose(fid);
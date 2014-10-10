function phi0 = rawread(fn)
fp = fopen(fn,'r');
nr = fscanf(fp,'%d',1);
nc = fscanf(fp,'%d',1);
phi0 = fscanf(fp,'%f',[nr nc])';
%phi0=fread(fp,nr*nc,'single=>single');
fclose(fp);


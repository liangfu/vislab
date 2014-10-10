function ksvd

clear all;clf;
winsize=8;


fp=fopen('../data/dict140x64.bin','r');
D=reshape(fread(fp,140*(winsize^2),'uint8'),winsize^2,[])';fclose(fp);





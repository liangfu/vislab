function [cy cx]=imcentroid(bw)
% IMCENTROID get image gravity center
[yy xx]=find(bw~=0);
cy=sum(yy)/size(yy,1);
cx=sum(xx)/size(xx,1);

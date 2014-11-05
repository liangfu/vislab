function neg_showhead()

channelLoc=[ 0  61 91;15  50 99;-16  51 104;0  42 113
             0 -61 91;15 -50 99;-16 -51 104;0 -42 113];
% channelLoc=[ -72  36 53;-72  24 67;-78  24 40;-82  12 56;
%              -72 -36 53;-75 -24 67;-78 -24 40;-82 -12 56];
values=[.1 .3 .7 1 -.1 -.3 -.7 -.9];

obj(1).verts=load('data/verts0.txt');
obj(1).faces=load('data/faces0.txt');
obj(1).norms=load('data/norms0.txt');
obj(1).color=-ones([size(obj(1).verts,1),1]);
%obj(1).color(1:122,1)=.8;

obj(2).verts=load('data/verts1.txt');
obj(2).faces=load('data/faces1.txt');
obj(2).norms=load('data/norms1.txt');
obj(2).color=-ones([size(obj(2).verts,1),1]);
obj(2).color=putcolor(obj,channelLoc(1:8,:),values(1:8));

labels=textread('data/labels.txt','%s');
locations=load('data/locations.txt');

channelLoc(:,3)=channelLoc(:,3)-60;
obj(1).verts(:,3)=obj(1).verts(:,3)-60;
obj(2).verts(:,3)=obj(2).verts(:,3)-60;
locations(:,3)=locations(:,3)-60;

figure(1);clf;

% % skull
% h1=patch('Vertices',obj(1).verts,'Faces',obj(1).faces,'FaceVertexCData',obj(1).color,'FaceColor','interp');
% set(h1,'FaceLighting','phong','EdgeLighting', 'phong', 'EdgeColor', 'none'); %, 'FaceColor', [.9 .9 .9]
% alpha(h1,.1)

hold on;

% % cortex
h2=patch('Vertices',obj(2).verts,'Faces',obj(2).faces,'FaceVertexCData',obj(2).color,'FaceColor','interp');
set(h2,'FaceLighting','phong','EdgeLighting','phong','EdgeColor','none'); % ,'FaceColor', [.9 .9 .9]
alpha(h2,.8)

caxis([-1 1]);
camzoom(1.2)
colormap(jet_modified);

plot3(locations(:,1),locations(:,2),locations(:,3),'.','markersize',16);
plot3(channelLoc(:,1),channelLoc(:,2),channelLoc(:,3),'r.','markersize',16);
text(locations(:,1),locations(:,2),locations(:,3),labels,'fontsize',10);

hold off;

light('position',[1 0 1]);
light('position',[-1 0 0]);
axis off;axis equal;%axis manual;
colorbar;

viewopt=2;
if viewopt==1,view(180,40);% side view
elseif viewopt==2,view(-90,90); % top view
elseif viewopt==3,view(-90,15); % back view
end
    
set(gcf,'toolbar','figure');

end

function color=putcolor(obj,loc,val)
verts=obj(2).verts;
faces=obj(2).faces;
color=obj(2).color;
R=minmax(loc');
extent=diff(R');
R(:,1)=R(:,1)-extent(:)*.2;
R(:,2)=R(:,2)+extent(:)*.2;

idx0=find(verts(:,1)>R(1,1));idx1=find(verts(:,1)<R(1,2));xidx=intersect(idx0,idx1);
idx0=find(verts(:,2)>R(2,1));idx1=find(verts(:,2)<R(2,2));yidx=intersect(idx0,idx1);
idx0=find(verts(:,3)>R(3,1));idx1=find(verts(:,3)<R(3,2));zidx=intersect(idx0,idx1);

idx=intersect(xidx,zidx);
idx=intersect(yidx,idx);

if mean(abs(loc(:,2)))>36,idx0=find(abs(verts(:,2))>36);idx=intersect(idx0,idx);end

if ~exist('scatteredInterp'),addpath('../utils');end
color(idx)=scatteredInterp(loc,val,verts(idx,:));

end

function c=jet_modified()
c=jet();c(1:2,:)=[1 1 1; 0 0 .5];
end

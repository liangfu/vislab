function writemodel_xml(model,filename)
    
[pathstr,basename,ext]=fileparts(filename);

PRECISION='%.8g ';

% clear filters
for ii=1:length(model.components)
comp=model.components{ii};
for jj=1:length(comp)
if length(comp(jj).filterid)>1
  ids=comp(jj).filterid(2:end);
  for kk=ids,model.filters(kk).w(:,:,:)=0;end
end
end
end

filtersw='';
for ii=1:length(model.filters)
  rows=size(model.filters(ii).w,1);
  cols=size(model.filters(ii).w,2)*size(model.filters(ii).w,3);
  data=model.filters(ii).w;
  data=reshape(permute(data,[3,2,1]),size(model.filters(ii).w,3),[]);
  filtersw=[filtersw '<_ type_id="opencv-matrix">\n'...
            sprintf('<rows>%d</rows>\n<cols>%d</cols>\n<dt>d</dt>\n<data>',rows,cols) ...
            sprintf(PRECISION,data) '</data></_>\n'];
end
biasw=zeros([length(model.bias),1]);
for ii=1:length(model.bias),biasw(ii)=model.bias(ii).w;end
anchors=zeros([length(model.defs),3]);
defs='';
for ii=1:length(model.defs),anchors(ii,:)=model.defs(ii).anchor;end
for ii=1:length(model.defs),defs=[defs '<_>' sprintf(PRECISION,model.defs(ii).w) '</_>\n'];end
indexers='';
for ii=1:length(model.components)
  comp=model.components{ii};
  component='';
  for jj=1:length(comp)
    if length(comp(jj).defid)>1,comp(jj).defid=comp(jj).defid(1);end
    component=sprintf([component '<part-%d>\n\t' ...
                       sprintf('<parentid>%d</parentid>\n\t',comp(jj).parent-1) ...
                       sprintf('<filterid>%s</filterid>\n\t',sprintf('%d ',comp(jj).filterid(1)-1)) ...
                       sprintf('<biasid>%s</biasid>\n\t',sprintf('%d ',comp(jj).biasid(1)-1)) ...
                       sprintf('<defid>%s</defid>\n',sprintf('%d ',comp(jj).defid-1)) ...
                       '</part-%d>\n'],jj-1,jj-1);
  end
  indexers=[indexers sprintf('<component-%d>%s</component-%d>',ii-1,component,ii-1)];
end
xmldata=sprintf(['<?xml version="1.0"?>\n' ...
                 '<opencv_storage>\n' ...
                 sprintf('<name>%s</name>\n',basename) ...
                 sprintf('<interval>%d</interval>\n',model.interval) ...
                 sprintf(['<thresh>' PRECISION '</thresh>\n'],model.thresh) ...
                 sprintf('<sbin>%d</sbin>\n',model.sbin) ...
                 '<norient>18</norient>\n' ...
                 '<flen>32</flen>\n'...
                 sprintf('<filtersw>%s</filtersw>\n',filtersw) ...
                 sprintf('<biasw>%s</biasw>\n',sprintf(PRECISION,biasw)) ...
                 sprintf('<anchors>%s</anchors>\n',sprintf('%d ',(anchors(:,1:2)-1)')) ...
                 sprintf('<defs>%s</defs>\n',defs) ...
                 sprintf('<indexers>%s</indexers></opencv_storage>\n',indexers)]);


fp=fopen(filename,'wt');
fprintf(fp,xmldata);
fclose(fp);



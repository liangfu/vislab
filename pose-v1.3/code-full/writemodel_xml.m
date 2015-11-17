function writemodel_xml(model,filename)
    
[pathstr,basename,ext]=fileparts(filename);

filtersw='';
for ii=1:length(model.filters)
  rows=size(model.filters(ii).w,1);
  cols=size(model.filters(ii).w,2)*size(model.filters(ii).w,3);
  filtersw=[filtersw '<_ type_id="opencv-matrix">\n'...
            sprintf('<rows>%d</rows>\n<cols>%d</cols>\n<dt>d</dt>\n<data>',rows,cols) ...
            sprintf('%f ',model.filters(ii).w) '</data></_>\n'];
end
biasw=zeros([length(model.bias),1]);
for ii=1:length(model.bias),biasw(ii)=model.bias.w;end
anchors=zeros([length(model.defs),3]);
defs='';
for ii=1:length(model.defs),anchors(ii,:)=model.defs(ii).anchor;end
for ii=1:length(model.defs),defs=[defs '<_>' sprintf('%f ',model.defs(ii).w) '</_>\n'];end
indexers='';
for ii=1:length(model.components)
  comp=model.components{ii};
  component='';
  for jj=1:length(comp)
    component=sprintf([component '<part-%d>\n\t' ...
                       sprintf('<parentid>%d</parentid>\n\t',comp(jj).parent) ...
                       sprintf('<filterid>%s</filterid>\n\t',sprintf('%d ',comp(jj).filterid)) ...
                       sprintf('<biasid>%s</biasid>\n\t',sprintf('%d ',comp(jj).biasid)) ...
                       sprintf('<defid>%s</defid>\n',sprintf('%d ',comp(jj).defid)) ...
                       '</part-%d>\n'],jj-1,jj-1);
  end
  indexers=[indexers sprintf('<component-%d>%s</component-%d>',ii-1,component,ii-1)];
end
xmldata=sprintf(['<?xml version="1.0"?>\n' ...
                 '<opencv_storage>\n' ...
                 sprintf('<name>%s</name>\n',basename) ...
                 sprintf('<interval>%d</interval>\n',model.interval) ...
                 sprintf('<thresh>%f</thresh>\n',model.thresh) ...
                 sprintf('<sbin>%d</sbin>\n',model.sbin) ...
                 '<norient>18</norient>\n' ...
                 '<flen>32</flen>\n'...
                 sprintf('<filtersw>%s</filtersw>\n',filtersw) ...
                 sprintf('<biasw>%s</biasw>\n',sprintf('%f ',biasw)) ...
                 sprintf('<anchors>%s</anchors>\n',sprintf('%d ',anchors)) ...
                 sprintf('<defs>%s</defs>\n',defs) ...
                 sprintf('<indexers>%s</indexers></opencv_storage>\n',indexers)]);


fp=fopen(filename,'wt');
fprintf(fp,xmldata);
fclose(fp);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% coefs = input_spline(rgbI,nknots)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function allows the user to fit a spline to image data.  It
% draws the input image rgbI, then asks the user to manually enter
% the closed curve that is to be fit by a B-spline.  It then has
% the user input where the knots of the B-spline are.  This should
% be done carefully, if there will be multiple B-splines
% representing the same object that must be compared.  The user
% should think of the points entered as correspondences.  The best
% (in the least-squares sense) B-spline is then fit to the entered
% curve with the assigned knot locations.  Finally, the user can
% move around the control points to fine tune the fit.  
%
% Input:
% rgbI: The image to display while drawing the spline.
% nknots: The number of curve pieces in the spline.
%
% Output: The control points of the spline input.  
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function varargout = input_spline(arg1,nknots)

global input_spline_coefs;
global input_spline_fig_h;
global input_spline_axis_h;
global input_spline_text_h;
global input_spline_button_h;
global input_spline_control_line_h;
global input_spline_control_points_h;
global input_spline_control_point_text_h;
global input_spline_knots_h;
global input_spline_control_line_h;
global input_spline_curve_h;
global input_spline_curve_p;
global selected_cp_h;
global input_spline_nknots;
global input_spline_image_h;

if strcmp(arg1,'finish'),
  set(gcbo,'UserData',1);
  set(input_spline_control_points_h,'ButtonDownFcn','');
  return;
elseif strcmp(arg1,'SelectControlPoint'),
  try
    select_control_point(gcbo);
  catch
    disp('error in select control point');
    set(input_spline_fig_h,'WindowButtonMotionFcn','','WindowButtonUpFcn','');
    set(input_spline_button_h,'UserData',1);
  end;
  return;
elseif strcmp(arg1,'MoveControlPoint'),
  try
    move_control_point(gcbo);
  catch
    disp('error in move control point');
    set(input_spline_fig_h,'WindowButtonMotionFcn','','WindowButtonUpFcn','');
    set(input_spline_button_h,'UserData',1);
  end;
  return;
elseif strcmp(arg1,'DeselectControlPoint'),
  try
    deselect_control_point(gcbo);
  catch
    disp('error in deselect control point');
    set(input_spline_fig_h,'WindowButtonMotionFcn','','WindowButtonUpFcn','');
    set(input_spline_button_h,'UserData',1);
  end
  return;
else
  rgbI = arg1;
end;

input_spline_nknots = nknots;
k = 4;

global L;
L = nknots;
setup_cubic_splines;

% Set up the figure window
[input_spline_fig_h,input_spline_axis_h,...
 input_spline_text_h,input_spline_button_h,input_spline_image_h] = ...
    set_up_window(rgbI);

old_db = get(input_spline_fig_h, 'DoubleBuffer');
state = uisuspend(input_spline_fig_h);
set(input_spline_fig_h,'DoubleBuffer','on');

% (1) Input a free form curve.  
[boundx,boundy] = input_free_form_curve(input_spline_text_h);

% Find points to approximate the spline

%s0 = 0:length(boundx);
%s = linspace(0,length(boundx),500);
%boundx = interp1(s0,[boundx;boundx(1)],s);
%boundy = interp1(s0,[boundy;boundy(1)],s);
[boundx,boundy] = reparameterize_spline(boundx,boundy,rgbI, ...
					input_spline_text_h, ...
					input_spline_axis_h);

% (2) Input what is a unit distance
ishappy = 0;
while ~ishappy,
  [u_x,u_y,u_i] = define_unit_distance(input_spline_text_h,boundx,boundy,nknots);
  
  % Find points on the curve evenly spaced between 
  [boundx,boundy] = sample_correspondence_uniformly(u_x,u_y,u_i,...
						    boundx,boundy,...
						    input_spline_text_h,nknots);
  
  % Fit a B-spline with nknots knots to the new [boundx,boundy]
  s = linspace(3,nknots+3,length(boundx)+1);
  s = s(1:length(boundx))';
  input_spline_coefs = fit_closed_b_spline(s,[boundx',boundy'],nknots)';
  
  cla;
  input_spline_image_h = image(rgbI);
  axis('image');
  hold on;
  
  % Draw the spline
  [input_spline_control_line_h,input_spline_curve_h,...
   input_spline_knots_h,...
   input_spline_curve_p,input_spline_control_points_h, ...
   input_spline_control_point_text_h]...
      = draw_annotated_spline(input_spline_coefs,nknots);
  
  
  set(input_spline_text_h,'string',{'Are you happy?','Left click to go on, right click to redo'});
  [tmp,tmp,b] = myginput(1);
  if b == 1,
    ishappy = 1;
  else,
    cla;
    input_spline_image_h = image(rgbI);
    axis('image');
    hold on;
    plot(boundx,boundy,'r','linewidth',2);
  end;
end  

% (3) Allow user to move around control points until he pushes the
% Done button
set(input_spline_button_h,'Callback','input_spline(''finish'')', ...
		  'Visible','on');
set(input_spline_text_h,'string',{'Improve the spline fit:',...
		    'Drag around the control points.',...
		    'Push the ''Done'' button when finished.\n'});
waitfor(input_spline_button_h,'UserData',1);
set(input_spline_text_h,'string','Done!');
set(input_spline_button_h,'Callback','', 'Visible','off');

set(input_spline_curve_h,'linewidth',3);
uirestore(state);
set(input_spline_fig_h, 'DoubleBuffer', old_db);

delete(input_spline_text_h);
delete(input_spline_button_h);

varargout = {input_spline_coefs};
clear global input_spline_coefs;
clear global input_spline_fig_h;
clear global input_spline_axis_h;
clear global input_spline_text_h;
clear global input_spline_button_h;
clear global input_spline_control_line_h;
clear global input_spline_control_points_h;
clear global input_spline_control_point_text_h;
clear global input_spline_knots_h;
clear global input_spline_control_line_h;
clear global input_spline_curve_h;
clear global input_spline_curve_p;
clear global selected_cp_h;
clear global input_spline_nknots;
clear global input_spline_image_h;

function [control_line_h,curve_h,knots_h,curve_p, ...
	  control_points_h,control_point_text_h]...
    = draw_annotated_spline(coefs,nknots)
hold on;
control_line_h = draw_control_line(coefs);
curve_p = spline_plot(coefs);
curve_h = plot(curve_p(:,1)',...
	       curve_p(:,2)',...
	       'r','HitTest','off');
knots_h = draw_knots(coefs);
[control_points_h,control_point_text_h] = draw_control_points(coefs);

function [control_points_h,control_point_text_h] = ...
    draw_control_points(coefs)

nknots = size(coefs,2)-3;

for i = 1:nknots,
  control_points_h(i) = plot(coefs(1,i+3),coefs(2,i+3),'gx');
  set(control_points_h(i),'UserData',i);
  control_point_text_h(i) = ...
      text(coefs(1,i+3),coefs(2,i+3)-3,...
	   num2str(i+2),'color','g','HitTest','off');
end;

set(control_points_h,'ButtonDownFcn','input_spline(''SelectControlPoint'')');
%set(control_points_h,'ButtonDownFcn','tmp_fcn');
function control_line_h = draw_control_line(coefs)

nknots = size(coefs,2)-3;

control_line_h = plot(coefs(1,[4:nknots+3,4]),...
		      coefs(2,[4:nknots+3,4]),...
		      'g-','HitTest','off');


function knots_h = draw_knots(coefs)
nknots = size(coefs,2)-3;
[x,y] = evaluate_spline_curve(repmat(coefs,[1,1,nknots]),3:nknots+2);
knots_h = plot(x,y,'yo','HitTest','off');

function [h_fig,h_axis,h_text,h_button,h_image] = set_up_window(rgbI)

% open/clear a figure
if ishandle(gcf), 
  clf; 
  h_fig = gcf;
else,
  figure;
  h_fig = gcf;
end;

% set figure color
set(h_fig,'color','k');

% Figure out the size of the screen
screen_size = get(0,'ScreenSize');
fig_position = [1+.1*screen_size(3),1+.1*screen_size(4),...
		.8*screen_size(3),.8*screen_size(4)];
set(h_fig,'Position',fig_position);

% Set up axes for drawing the image and plotting the spline
h_axis = axes('position',[.02,.3,.96,.68]);

% Draw the current image
h_image = image(rgbI); axis('image');

plot_position = get(h_axis,'position');

a = 1 - plot_position(4) - plot_position(2);
b = plot_position(4);
c = 1/3 * (1 - 4 * a - b);

% Set up a text box for entering directions
text_position = [plot_position(1),2*a+c,...
		 plot_position(3),2*c];
h_text = uicontrol('style','text',...
		   'Units','normalized',...
		   'position',text_position,...
		   'BackgroundColor',[51,51,51]/255,...
		   'ForegroundColor',[255,255,255]/255,...
		   'String','Directions',...
		   'FontName','Times',...
		   'FontSize',14,...
		   'HorizontalAlignment','center');

% set up a button for indicating when done
button_position = [((1 - plot_position(1)*2)*fig_position(3) - ...
		    60)/2+plot_position(1)*fig_position(3),
		   a * fig_position(4)+ (c * fig_position(4) - 20)/2,
		   60; 20];
%button_position = [plot_position(1),a,4*c,c];
h_button = uicontrol('style','pushbutton',...
		      'Units','pixels',...
		      'Position',button_position,...
		      'String','Done',...
		      'Callback','',...
		      'FontName','Times',...
		      'FontSize',14,...
		      'HorizontalAlignment','center',...
		      'Visible','off',...
		      'UserData',0);

function [x,y] = reparameterize_spline(x0,y0,rgbI,h_text,h_axis)

global input_spline_image_h;

set(h_text,'string','Please wait ... Reparameterizing the spline');
drawnow;

s0 = 0:length(x0);
s = linspace(0,length(x0),500);
x = interp1(s0,[x0;x0(1)],s);
y = interp1(s0,[y0;y0(1)],s);

% Draw the boundary
axes(h_axis);
cla;
input_spline_image_h = image(rgbI);
axis('image');
is_hold = ishold;
hold on;
plot(x,y,'r.');
if is_hold==0,
  hold off;
end;

function [u_x,u_y,u_i] = define_unit_distance(h_text,boundx,boundy,nknots)

set(h_text,'string',{'Define unit distance along the B-spline curve:',...
		    'Left click to enter a point unit distance from the previous point (this point will be a knot of the spline).',...
		    'The first knot is the first point clicked when free drawing the curve.',...
		    sprintf('There should be a total of %d points entered.',nknots)});

% First knot is the first point clicked
u_x(1) = boundx(1); u_y(1) = boundy(1); u_i(1) = 1;
is_hold = ishold;
hold on;
plot(u_x(1),u_y(1),'mo');
text(u_x(1),u_y(1)-3,'1','color','m');

for i = 2:nknots,
  [tmpx,tmpy] = myginput(1);
  
  % find the closest point on the contour
  [mindist,minind] = min((tmpx-boundx).^2 + (tmpy-boundy).^2);
  u_x(i) = boundx(minind);   u_y(i) = boundy(minind);
  u_i(i) = minind;
  
  % plot the point
  plot(u_x(i),u_y(i),'mo');
  text(u_x(i),u_y(i)-3,num2str(i),'color','m');
  drawnow;
  
end;

u_x(end+1) = boundx(end); u_y(end+1) = boundy(end);
u_i(end+1) = length(boundx);

set(h_text,'string','Define unit distance along the B-spline curve: DONE!');
if is_hold == 0,
  hold off;
end;

function [new_boundx,new_boundy] = ...
    sample_correspondence_uniformly(u_x,u_y,u_i,boundx,boundy,h_text,nknots)

% number of samples between each correspondence is the maximum
% number of boundary points between any two correspondences
interval_length = diff(u_i)+1;
[max_length,minind] = max(interval_length);
nsamples = max(15,max_length);

set(h_text,'string','Interpolating uniformly between the input knots.');

for i = 1:nknots,
  
  start_i = u_i(i);
  end_i = u_i(i+1);
  l = end_i - start_i + 1;
  s = start_i:end_i;
  si = linspace(start_i,end_i-l/nsamples,nsamples);
  px = boundx(start_i:end_i);
  py = boundy(start_i:end_i);
  new_boundx((i-1)*nsamples+1:i*nsamples) = interp1(s,px,si);
  new_boundy((i-1)*nsamples+1:i*nsamples) = interp1(s,py,si);

end;

function select_control_point(h)

global selected_cp_h;
global selected_cp_ind;
global input_spline_fig_h;
global input_spline_control_points_h;

selected_cp_h = h;
selected_cp_ind = find(selected_cp_h == input_spline_control_points_h);
set(input_spline_fig_h,'WindowButtonMotionFcn',...
		  'input_spline(''MoveControlPoint'')',...
                  'WindowButtonUpFcn',...
		  'input_spline(''DeselectControlPoint'')');

function deselect_control_point(h)

global selected_cp_h;
global selected_cp_ind;
global input_spline_fig_h;
global input_spline_coefs;
global input_spline_knots_h;
global input_spline_control_point_text_h;
global input_spline_nknots;

% redraw the knot locations
[x,y] = evaluate_spline_curve(repmat(input_spline_coefs,[1,1, ...
		    input_spline_nknots]),3:input_spline_nknots+2);
set(input_spline_knots_h,'XData',x,'YData',y);

% Move the text for the current cp
x = get(selected_cp_h,'XData');
y = get(selected_cp_h,'YData');
set(input_spline_control_point_text_h(selected_cp_ind),'Position',[x,y-3,0]);

set(input_spline_fig_h,'WindowButtonMotionFcn','','WindowButtonUpFcn','');

selected_cp_h = [];

function move_control_point(h)
global selected_cp_h;
global selected_cp_ind;
global input_spline_coefs;
global input_spline_fig_h;
global input_spline_axis_h;
global input_spline_image_h;
global input_spline_curve_p;
global input_spline_curve_h;
global input_spline_control_line_h;
global input_spline_control_points_h;
global input_spline_nknots;

k = 4;

[imageHandle,x,y] = over_image(input_spline_axis_h);
if imageHandle ~= 0,
  if imageHandle == input_spline_image_h;

    % set control point to the current mouse location
    new_cp = [x;y];
    set(selected_cp_h,'XData',new_cp(1),'YData',new_cp(2));
    input_spline_coefs(:,selected_cp_ind+3) = new_cp;
    
    % end conditions: some of the control points are redundant, so
    % replicate changes in redundant control points
    if input_spline_nknots-selected_cp_ind < 3,
      input_spline_coefs(:,-input_spline_nknots+selected_cp_ind+3) = ...
	  new_cp;
    end;
    
    % redraw the lines connecting the control points
    set(input_spline_control_line_h,'XData',...
		      input_spline_coefs(1,[4:input_spline_nknots+3,4]),...
		      'YData',input_spline_coefs(2,[4:input_spline_nknots+3,4]));
    
    % redraw the spline
    input_spline_curve_p = spline_plot(input_spline_coefs);
    set(input_spline_curve_h,'XData',input_spline_curve_p(:,1)',...
		      'YData',input_spline_curve_p(:,2)');
  end;
end;

function [imageHandle,x,y] = over_image(axesHandle)

% Return the index of which image we are over, and return a 0 if we
% aren't above an image.

imagefound = findobj(axesHandle, 'type', 'image');
if isempty(imagefound)
   imageHandle=0; x=0; y=0;
   return
end
% Make sure that the Image's Button Down & Up functions will queue
set(imagefound, 'Interruptible', 'off', 'BusyAction', 'Queue');
axHandle = get(imagefound, 'Parent');
axPosition = get(axHandle, 'Position');
axCurPt = get(axHandle, 'CurrentPoint');

% See if we are above the desired axes
imageHandle = 0;  
XLim = get(axHandle, 'XLim');
YLim = get(axHandle, 'YLim');
pt = axCurPt;
x = pt(1,1); y = pt(1,2);
if x>=XLim(1) & x<=XLim(2) & y>=YLim(1) & y<=YLim(2)
  imageHandle = imagefound;
else,
  x = 0; y = 0;
end;


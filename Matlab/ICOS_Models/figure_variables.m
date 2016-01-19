function figure_variables
% Draws a figure identifying the variables on a mirror surface
B = 0.4;
n = 7;
phi = pi/n;
rmin = B/(2*sin(phi));
r = 2*rmin;
w = r*sin(phi);
h = r*cos(phi);
clf;
draw_circle([ w, h], B/2, 'g');
draw_circle([-w, h], B/2, 'g');
H1 = [
  draw_radius(-phi, r, r/2, [-B/4 0], '\it{r}');
  draw_radius(phi, r, r/2, [B/4 0], '\it{r}','HorizontalAlignment','Right');
  draw_radius(0, h, r*0.65, [B/10 0], '\it{h}') ];
pp = w/8;
plot([-w,w],[h h],'k',[-pp -pp 0], [h, h-pp, h-pp],'k');
H2 = text(w/2, h+0.6*B, '\it{w}', 'HorizontalAlignment','Center');
%draw_arrow(ax, O, D, headlen, limrange)
draw_arrow(gca, [-B/2,h+0.6*B],[0,h+0.6*B], B/8, [-0.55*B,B/4]);
draw_arrow(gca, [w+B/2,h+0.6*B],[w,h+0.6*B], B/8, [-0.55*B,B/4]);

H3 = text(-w-B,h,'\it{B}','HorizontalAlignment','Center');
draw_arrow(gca, [-w-B,h+B], [-w-B,h+B/2], B/8, [-B,B/4]);
draw_arrow(gca, [-w-B,h-B], [-w-B,h-B/2], B/8, [-B/4,B]);

H4 = text(B/8, 0, '\it{O}');


H5 = text(-(r*0.6)*sin(phi/2), (r*0.6)*cos(phi/2), '\phi', ...
  'HorizontalAlignment','Center');
draw_arc_arrow(gca, [0 0], r/2, -[phi/2 phi], B/8);
draw_arc_arrow(gca, [0 0], r/2, -[phi/2 0], B/8);

set([H1;H2;H3;H4;H5],'FontUnits','Normalized','FontSize',1/12);
set(gca,'DataAspectRatio',[1 1 1],'XTick',[],'YTick',[]);
ylim([-B/2 h+3*B/2]);
xlim((w+3*B/2)*[-1 1]);
shg;
print_figure(gca,'figure_variables', [3 3], [H1;H2;H3;H4;H5]);

function draw_circle(O, r, c)
sth = linspace(0,2*pi,21);
x = O(1)+r*sin(sth);
y = O(2)+r*cos(sth);
plot(x,y,c);
hold on;

function h = draw_radius(phi, r, lblr, offset, varargin)
X = [0 r*sin(phi)];
Y = [0 r*cos(phi)];
hold on;
plot(X,Y,'k');
if nargin >= 5
  h = text(lblr*sin(phi)+offset(1), lblr*cos(phi)+offset(2), ...
    varargin{:});
  % set(h,'FontUnits','Normalized');
end

function print_figure(ax, fname, target_dim, t)
  target_width=target_dim(1);
  target_height=target_dim(2);
  
  target_aspect_ratio = target_height/target_width;
  set(ax,'units','inches');
  P = get(ax,'position');
  source_aspect_ratio = P(4)/P(3);
  if source_aspect_ratio >= target_aspect_ratio
    % fit to target height
    height = target_height;
    width = height/source_aspect_ratio;
    print_scale = height/P(4);
  else
    width = target_width;
    height = width * source_aspect_ratio;
    print_scale = width/P(3);
  end
  
  set(gcf,'InvertHardcopy','on');
  set(gcf,'PaperUnits', 'inches');
  papersize = get(gcf, 'PaperSize');
  % Center the image on the page:
  left = (papersize(1)- width)/2;
  bottom = (papersize(2)- height)/2;
  myfiguresize = [left, bottom, width, height];
  set(gcf,'PaperPosition', myfiguresize);
  
  % Now change the font size to map the size we will be scaled to
  tpixels = zeros(size(t));
  for i=1:length(t)
    set(t(i),'fontunits','pixels');
    tpixels(i) = get(t(i),'fontsize');
    set(t(i),'fontsize', tpixels(i)*print_scale);
  end
  
%   h12w = get(h12,'linewidth');
%   set(h12,'linewidth',h12w*print_scale);
%   
%   h23w = get(h23,'linewidth');
%   set(h23,'linewidth',h23w*print_scale);
  
  set(ax,'xlimmode','manual','ylimmode','manual','zlimmode','manual');
  
  print(fname,'-dpng','-r300');
  
  for i=1:length(t)
    set(t(i),'fontsize', tpixels(i));
  end
  % set(h12,'linewidth',h12w);
  % set(h23,'linewidth',h23w);

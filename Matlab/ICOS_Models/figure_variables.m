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

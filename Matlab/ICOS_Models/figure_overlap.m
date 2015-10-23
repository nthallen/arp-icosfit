function figure_overlap
% figure_overlap
% Draws a figure to explain the overlap calculation
n = 13;
B = .4;
w = B/2;
phi = pi/n;
rmin = w/sin(phi);
rmax = 2*rmin;
N = 1:n;
phis = 2*N*phi; % angle of spot centers, starting at vertical
O = rmin*[sin(phis') cos(phis')];
for i=1:n
  draw_circle(O(i,:), w, 'g');
end

phis3 = 2*N*phi; % same for outer spots
O3 = rmax*[sin(phis3') cos(phis3')];
for i=1:n
  draw_circle(O3(i,:), w, 'b');
end

% Indicate spot width, using last spot
phi4 = 2*N(end)*phi;
O4 = rmin*[sin(phi4) cos(phi4)];
dy = B;
draw_arrow(gca,O4+[-B B], O4+[-B/2 B],B/4,[-B B/4]);
draw_arrow(gca,O4+[B B], O4+[B/2 B],B/4,[-B B/4]);
hold on;
h = text(O4(1),O4(2)+B,'\it{B}','HorizontalAlignment','Center');
set(h,'FontUnits','Normalized');

draw_radius(2*7*phi, rmin, rmin/2, [B/4 0], '\it{r_{min}}');
draw_radius(2*8*phi, rmax, .75*rmax, [-B/3,B/4],'\it{r_{max}}', ...
  'horizontalalignment','center');

% Draw angle 2phi around the 3rd spot
phis2 = (2*3 + [-1 1])*phi;
E = (rmax+w)*[sin(phis2') cos(phis2')];
plot([0 E(1,1)],[0 E(1,2)],'k',[0 E(2,1)],[0 E(2,2)],'k');
draw_arc_arrow(gca, [0 0], 1.5*rmin, phis2(1)+[-phi,0], B/4);
draw_arc_arrow(gca, [0 0], 1.5*rmin, phis2(2)+[phi,0], B/4);

h = text(1.5*rmin*sin(phis(3)), 1.5*rmin*cos(phis(3)), '2\phi');
set(h,'horizontalalignment','center','FontUnits','Normalized');

hold off;
set(gca,'DataAspectRatio',[1 1 1],'XTick',[],'YTick',[]);
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
  set(h,'FontUnits','Normalized');
end

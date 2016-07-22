function draw_targets(P,target_spacing)
% draw_targets(P)
% draw_targets(P,target_spacing)
%   target_spacing is a 2x1 specifying distance from the
%   back face of the first ICOS mirror in cm.
%   target_spacing(1) is the position of the first target
%   target_spacing(2) is the position of the target near the
%   back surface of the first ICOS mirror.
% Draws two targets, one for the Herriott mirror and one
% for the first ICOS mirror.
scale = 1/2.54;
if nargin < 2
  target_spacing = [P.herriott_spacing,0];
end
r0 = .472*2.54; % slot minimum radius (from optical axis)
rs = .098*2.54; % slot radius
rs = P.beam_diameter/2;
if ~isfield(P,'r1')
  error('ICOS_Model5 is not supported');
  target_radius = P.r; % ICOS_Model5
end
target_radiusI = P.r1; % ICOS_Model6
target_radiusR = P.Hr;
m = P.injection_scale;
clf;
D = [1,-P.dy*m,P.dz*m];
IP0 = [-P.CT1, m*P.y0, m*P.z0];
IP = IP0 - target_spacing(2)*D;
HP = IP0 - target_spacing(1)*D;
rot = atan2d(HP(3),HP(2));
M = [1 0 0
  0 cosd(rot) -sind(rot)
  0 sind(rot)  cosd(rot)];
IPr = IP*M;
HPr = HP*M;
txtI = sprintf('Position %.2f cm', target_spacing(2));
txtR = sprintf('Position %.2f cm', target_spacing(1));

Ioffset = target_radiusI + target_radiusR + 0.5;
% Ioffset = 0;
draw_target('ICOS', txtI, IPr, target_radiusI, Ioffset, 0, scale, rs);
draw_target('RIM', txtR, HPr, target_radiusR, 0, rs, scale, rs);
draw_slot(r0, rs, target_radiusR, scale);
hold off;
shg;
set(gca,'xtick',[],'ytick',[]);
drawnow;
yl = ylim(gca);
xl = xlim(gca);
ylim(gca,yl);
xlim(gca,xl);
y_x = (yl(2)-yl(1))/(xl(2)-xl(1));
fp = get(gcf,'position');
y_top = fp(2)+fp(4);
fp(4) = fp(3)*y_x;
% fp(2) = fp(2)-fp4+fp(4);
% fp(4) = fp4;
% fp
set(gcf,'position',fp);
drawnow;
fp2 = get(gcf,'position');
if fp2(4) < fp(4)
  fp2(3) = fp2(4)/y_x;
end
fp2(2) = y_top-fp2(4);
set(gcf,'position',fp2);
set(gca,'position',[0 0 1 1 ]);
pp = get(gcf,'paperposition');
pp(3) = xl(2)-xl(1);
pp(4) = yl(2)-yl(1);
set(gcf,'paperposition',pp);
drawnow;
h = findobj(gca,'type','text');
set(h,'fontunits','normalized');

function draw_target(ttl, subttl, tgt, r, zoffset, rs, s, beam_r)
% draw_target(ttl, tgt, r, zoffset)
if nargin < 8
  beam_r = 0.2;
end
res = 360;
t0 = atan(rs/r);
theta = pi+linspace(t0,2*pi-t0,res);
x = r*cos(theta);
y = r*sin(theta)+zoffset;
ax = [-r-0.3,r+0.3];
z = [0 0];
ch = [-0.2,0.2];
plot(x*s,y*s,'k', ...
  -(tgt(2)+ch)*s, (tgt(3)+zoffset+z)*s, 'k', ...
  -(tgt(2)+z)*s, (tgt(3)+zoffset+ch)*s, 'k', ...
  ax*s,(zoffset+z)*s, 'k', ...
  z*s, (zoffset + ax)*s, 'k');
hold on;
res = 20;
theta = linspace(t0,2*pi-t0,res);
x = (-tgt(2)+beam_r*cos(theta))*s;
y = (tgt(3)+beam_r*sin(theta)+zoffset)*s;
plot(x,y,'k');
x = (beam_r*cos(theta))*s;
y = (beam_r*sin(theta)+zoffset)*s;
plot(x,y,'k');
set(gca,'DataAspectRatio',[1 1 1],'visible','off');
h = [ text(-0.2*s, (0.4+zoffset)*s, ttl);
  text(-0.2*s, (-0.4+zoffset)*s, subttl) ];
set(h,'horizontalalignment', 'right');
set(h(2),'verticalalignment', 'top');

function draw_slot(r0, rs, rmax, s)
res = 180;
theta = linspace(0,pi,res);
x0 = sqrt(rmax^2-rs^2);
x = -[x0, r0-rs*sin(theta), x0];
y = [rs, rs*cos(theta), -rs];
plot(x*s,y*s,'k');

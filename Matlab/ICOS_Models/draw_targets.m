function draw_targets(P,scale)
% draw_targets(P)
% draws two targets, one for the Herriott mirror and one
% for the first ICOS mirror.
if nargin < 2
  scale = 1;
end
r0 = .472*2.54;
rs = .098*2.54;
m = P.injection_scale;
clf;
HP = [-P.herriott_spacing, m*P.y0, 0];
D = [1,-P.dy*m,P.dz*m];
IP = HP+P.herriott_spacing*D;
txt = sprintf('Herriott Spacing %.2f cm', P.herriott_spacing);

Ioffset = P.r * 2 + 0.5;
% Ioffset = 0;
draw_target('ICOS', txt, IP, P.r, Ioffset, 0, scale);
draw_target('RIM', txt, HP, P.r, 0, rs, scale);
draw_slot(r0, rs, P.r, scale);
hold off;
shg;
set(gca,'xtick',[],'ytick',[]);
yl = ylim(gca);
xl = xlim(gca);
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

function draw_target(ttl, subttl, tgt, r, zoffset, rs, s)
% draw_target(ttl, tgt, r, zoffset)
res = 360;
t0 = atan(rs/r);
theta = linspace(t0,2*pi-t0,res);
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
x = (-tgt(2)+0.2*cos(theta))*s;
y = (tgt(3)+0.2*sin(theta)+zoffset)*s;
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
x = [x0, r0-rs*sin(theta), x0];
y = [rs, rs*cos(theta), -rs];
plot(x*s,y*s,'k');

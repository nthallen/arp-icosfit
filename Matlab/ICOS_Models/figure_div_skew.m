%%
f = figure;
set(f,'menubar','none','position',[107 80 1091 711]);
%%
xlen = 2;
x0 = 0.5;
y0 = 0.3;
dx = 1;
dy = 0.1;
dz = 0.1;
dT = 0.1;
FS = 16;
%%
clf
ax = axes;
view(3);
xlim(ax,[0 xlen]);
ylim(ax,2*dT*[-1,1]);
zlim(ax,[-dT y0+dy+dT]);
set(ax,'DataAspectRatio',[1 1 1],'visible','off');
set(f,'color',[1 1 1]);
set(ax,'units','pixels');
drawnow;
%
xpos = get(ax,'position');
set(f,'units','pixels');
fpos = get(f,'position');
menuheight=80;
if xpos(4)/xpos(3) > (fpos(4)-menuheight)/fpos(3)
  nfpos = [fpos(1),fpos(2),(fpos(4)-menuheight)*xpos(3)/xpos(4),fpos(4)];
  % fprintf(1,'axes taller than figure\n');
else
  nfpos = [fpos(1),fpos(2),fpos(3),fpos(3)*xpos(4)/xpos(3)+menuheight];
  % fprintf(1,'axes wider than figure\n');
end
set(f,'position',nfpos);
% fprintf(1, ' fpos = [ %s]\n', sprintf('%d ',floor(fpos)));
% fprintf(1, 'nfpos = [ %s]\n', sprintf('%d ',floor(nfpos)));
nxpos = [0 0 nfpos(3) nfpos(3)*xpos(4)/xpos(3)];
% fprintf(1, ' xpos = [ %s]\n', sprintf('%d ',floor(xpos)));
% fprintf(1, 'nxpos = [ %s]\n', sprintf('%d ',floor(nxpos)));
% fprintf(1, 'aspects = [ %s]\n', sprintf('%.1f ', ...
%   [fpos(3)/(fpos(4)-menuheight) nfpos(3)/(nfpos(4)-menuheight) xpos(3)/xpos(4) nxpos(3)/nxpos(4)]));
set(ax,'position',nxpos);
%
draw_arrow3(ax, [0,0,0],[xlen,0,0],0.1);
text(xlen/2,0,-dT,'\it{x}','parent',ax,'fontsize',FS);
%%
H = draw_arrow3(ax, [x0,0,y0],[x0+dx,dz,y0+dy]);
set(H,'color',[1 0 0]);
%%
draw_arrow3(ax, [x0,0,0],[x0,0,y0],0.1);
text(x0-dT,0,y0/2,'\it{r}','parent',ax,'fontsize',FS);
%%
draw_arrow3(ax, [x0,0,y0],[x0+dx,0,y0]);
draw_arrow3(ax, [x0+dx,0,y0],[x0+dx,0,y0+dy],dy/3);
draw_arrow3(ax, [x0+dx,0,y0],[x0+dx,dz,y0],dz/3);
draw_arrow3(ax, [x0+dx,0,y0],[x0+dx,dz,y0+dy],dz/3);
hold(ax,'on');
plot3(ax,(x0+dx)*[1 1 1],[dy dy 0],[y0 y0+dy y0+dy],'c');
%%
text(x0+dx/2,0,y0-dT,'dx=1','parent',ax,'fontsize',FS,'horizontalalignment','center');
%%
text(x0+dx+dT/2,0,y0+dy/3,'\it{d}','parent',ax,'fontsize',FS);
%%
text(x0+dx-dT/2,dz,y0,'\it{s}','parent',ax,'fontsize',FS);
%%
text(x0+dx/2,0,y0+dT/2,'\theta','parent',ax,'fontsize',FS);
%%
text(xlen/2+dT*3,0,-dT,'$$\mathrm{tan}\theta=\sqrt{d^2+s^2}$$','parent',ax,'fontsize',FS,'interpreter','latex');
%%
print(f, '-dpng','-r300','figure_div_skew.png');
%%
xlen = 2;
dT = 0.1;
FS = 16;
f = figure;
% set(f,'menubar','none','position',[107 80 1091 711]);
%%
clf
ax = axes;
view(3);
%xlim(ax,[0 xlen]);
%ylim(ax,2*dT*[-1,1]);
%zlim(ax,[-dT y0+dy+dT]);
set(ax,'DataAspectRatio',[1 1 1]);
set(f,'color',[1 1 1]);
set(ax,'units','pixels');
drawnow;
%
draw_arrow3(ax, [0,0,0],[xlen,0,0],0.1);
text(xlen/2,0,-dT,'\it{x}','parent',ax,'fontsize',FS);
xlim([-.2 xlen+1]);
ylim([-1 1]);
zlim([-.7 .5]);
% Now pick a good ray. I'd like the closest approach to be at
% xlen/2. Let's call that distance rf (for radius at focus)
% According to ICOS Theory doc, x = -rd/(s^2+d^2), and we want
% x = xlen/2. The radius at x will be rf = sr/sqrt(s^2+d^2)
r0 = 0.4;
rf = 0.1;
xf = xlen/2;
d0 = (rf^2-r0^2)/(r0*xf);
s0 = sqrt(-d0*rf^2/(r0*xf));
% 
x0 = 0;
y0 = 0;
z0 = r0;
dx = 1;
dy = s0;
dz = d0;
x = 0:.1:xlen;
y = y0 + x*dy;
z = z0 + x*dz;
hold(ax,'on');
plot3(ax,x,y,z,'r');
%
r = sqrt(y.^2+z.^2);
d = (y*dy+z*dz)./r;
s = (y*dz - z*dy)./r;
%%
view([-90 0]);
%%
view(3);
%%
shg;
arrowhead = 0.1;
i = 4;
xvec = [1 0 0];
ldir = [1 dy dz];
for i=1:length(x)
  rpos = [x(i) y(i) z(i)];
  rvec = [0 y(i) z(i)]/sqrt(y(i)^2+z(i)^2);
  svec = cross(xvec,rvec);
  xpos = [x(i) 0 0];
  di = sum(ldir.*rvec);
  si = sum(ldir.*svec);
  h1 = draw_arrow3(ax, xpos, rpos,arrowhead);
  set(h1,'color',[0 0 1]);
  fprintf(1,'rvec: %f svec: %f  xvec: %f  di: %f   si: %f\n', ...
    sqrt(sum(rvec.^2)), sqrt(sum(svec.^2)), sqrt(sum(xvec.^2)), di, si);
  %% divergence vector
  h2 = draw_arrow3(ax, rpos, rpos+di*rvec,arrowhead);
  set(h2,'color',[1 1 0]);
  %% skew vector
  h3 = draw_arrow3(ax, rpos, rpos+si*svec, arrowhead);
  set(h3,'color',[1 0 1]);
  %%
  h4 = draw_arrow3(ax, rpos, rpos+xvec, arrowhead);
  legend('x','laser','r','d','s');
  pause;
  delete([h1 h2 h3 h4]);
end

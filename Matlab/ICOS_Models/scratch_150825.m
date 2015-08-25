%%
d = -.2;
s = .02;
r = 3;
f = -5;

x = (0:.1:18)';
rx = sqrt((r+d*x).^2 + (x*s).^2);
dx = (r*d+x*(d^2+s^2))./rx;
figure; plot(x,rx);
%%
plot(x,dx); shg
%%
sx = s*r./rx;
plot(x,sx); shg
%%
th = atand(sqrt(dx.^2+sx.^2));
plot(x,th); shg
%%
f = -15;
d1x = dx - rx/f;
plot(x,dx,x,d1x); shg
%%
% d(dx)/dx
drx = (d*(r+d*x)+s^2*x)./rx;
figure; plot(x,rx,x,drx);
%%
ddx = (rx*(d^2+s^2)-(r*d+x*(d^2+s^2)).*drx)./(rx.^2);
figure; plot(x,dx,x,ddx);
legend('dx','ddx');
grid;
%%
dd1x = ddx - drx/f;
figure; plot(x,dx,x,ddx,x,d1x,x,dd1x,x,drx); grid;
legend('dx','ddx','d1x','dd1x','drx');
%%
% Check ddx:
x1 = x(1:end-1)+diff(x)/2;
ddxd = diff(dx)./diff(x);
figure; plot(x,ddx,x1,ddxd); legend('ddx','ddxd');
%%
V = [
  -(d^2+s^2)^2/f
  -3*r*d*(d^2+s^2)/f
  -r^2*(3*d^2+s^2)/f
  r^2*(s^2-r*d/f)
  ];
dd1xf = polyval(V,x);
figure; plot(x,d1x,x,dd1xf); legend('d1x','dd1xf');
grid;
%%
rts = roots(V);
rts = rts(imag(rts)==0);
rts = rts(rts > 0);
rts = min(rts); % rts should be where the radius == rmax

%%
% Testing equations for R1==R2 portion of the paper
L = 35;
r = 1;
W = .4;
rd = .1;
th = 15;
C = 2000;
%%
% R as a function of r over 0:3.8 cm (3" diameter)
% (It makes more sense to treat r as a function of R, which we
%  do in later cases, but this works). The r we are determining
% here is the upper limit for focusing.
r0 = sqrt(L*rd*tand(th));
r = r0:0.01:3.8;
r4 = r.^4;
r2 = r.^2;
Lrt = L*rd^2*tand(th)^2;
r4mL2rt = sqrt(max(0,r4-L*Lrt));
R1 = (r4 - r2.*r4mL2rt)/Lrt;
R2 = (r4 + r2.*r4mL2rt)/Lrt;
%
% plot(R1,r,R2,r);
% % xlim([0 max(r)]);
% xlabel('r');
% ylabel('R');
% shg;
%%
% Focusing Criteria
logX = false;
if logX
  drawplot=@semilogx;
  Xmin = L/2.25;
else
  drawplot=@plot;
  Xmin = 0;
end
Rup = [flip(R1) R2(2:end)];
rup = [flip(r) r(2:end) ];
r2k = interp1(Rup,rup,2000);
n2k = find(Rup<2000,1,'last');
n2100 = find(Rup<2100,1,'last');
Rbound = [ L/2    Rup(1:n2k) 2000 2000 L/2 L/2 ];
rbound = [ max(r) rup(1:n2k) r2k  0    0   max(r) ];
drawplot(Rup(1:n2100), rup(1:n2100));
hold on;
drawplot([2000 2000], [0 max(r)],'k',[L/2 L/2], [0 max(r)],'k');
% Fill in the focus region
h = fill(Rbound,rbound,'c');
xlabel('R, cm');
ylabel('r, cm');
hold off;
xlim([Xmin 2150]);
title(sprintf('Focusing Criteria, L=%d',L));
shg;
%%
% Add overlap criteria, Interleave = 1
n = ceil(C/(2*L));
phi = pi/n;
ro = W/(2*sin(phi));
delete(h);
hold on;
if ro < r0
  % one region
  Rbound = [ L/2    Rup(1:n2k) 2000 2000 L/2 L/2 ];
  rbound = [ max(r) rup(1:n2k) r2k  ro   ro  max(r) ];
  h = fill(Rbound,rbound,'c');
else
  % two regions
  Ro1 = interp1(r,R1,ro);
  nro1 = find(r>ro,1);
  Rbound1 = [ Ro1 R1(nro1:end) L/2    L/2    Ro1 ];
  rbound1 = [ ro  r(nro1:end)  max(r) ro     ro ];
  Ro2 = interp1(r,R2,ro);
  r2000 = interp1(R2,r,2000);
  nro2 = find(r>ro & R2 < 2000);
  Rbound2 = [Ro2 R2(nro2) 2000  2000 Ro2];
  rbound2 = [ro  r(nro2)  r2000 ro   ro ];
  h = [fill(Rbound1, rbound1, 'c') fill(Rbound2, rbound2, 'c')];
end
drawplot([Xmin Rup(n2100)], [ro ro], 'r');
hold off;
title(sprintf('Focusing and Overlap Criteria, L=%d',L));
shg;
%%
% Calculate R limits based on overlap criteria: not enough that r
% is large enough: w must also be small enough, so w/r < sin(pi/n)
Ro3 = L/(1+cos(pi/n));
Ro4 = L/(1-cos(pi/n));
%%
% Calculate RIM length based on interleave 1:
LR = W^2/(4*rd*tand(th)*sin(phi));
fprintf(1, 'RIM length: %.1f cm\n', LR);
%% Add 75 cm line
hold on;
plot([75 75], [0 max(r)], 'g');
hold off;
shg;
%%
figure;
%%
% For a given L, r(L) (chosen as the smallest r that will meet the
% overlap criteria with w = W/2), these graphs show:
%   a curve showing skew*r for vs R
%   a blue box indicating values of skew*r that both
%     -are below the limit for focusing on (2 mm)^2 detector, and
%     -are above the limit for meeting the overlap criteria
%   The former is really irrelevant
for L=50:10:100
  clf;
  r_min = W/(2*sin(2*L*pi/C));
  r = r_min;
  g = 0.01;
  maxR = (g^2+(r+W).^2)/(2*g);

  % maxR = 350;
  R = (L/2):.1:500;
  s = sqrt(2*R/L-1)*r./R;

  plot(R,s*r); xlabel('Radius of Curvature, cm'); ylabel('skew \times r');
  title(sprintf('L=%.0f, r=%.1f',L,r));
  %
  Wcon = W*r/(2*L);
  Dcon = rd*tand(th);
  if Wcon < Dcon
    color = 'c';
  else
    color = 'r';
  end
  x = [0, maxR, maxR, 0, 0];
  y = [Wcon, Wcon, Dcon, Dcon, Wcon];
  hold on;
  h = fill(x,y,color);
  alpha(h,.25);
  hold off;
  shg;
  pause;
end
%%
% Let me try to look at overlap criteria more directly
% w > W/2
% sin(2*L*pi/C) > w/r = (L/R)*sqrt(2*R/L-1)
for L=50:10:90
  clf;
  R = (L/2):.1:2000;
  sin2Lpi_C = sin(2*L*pi/C)*ones(size(R));
  rr = sqrt(2*R/L-1)*L./R;
  plot(R,rr,R,sin2Lpi_C);
  xlabel('R'); ylabel('w/r');
  legend('w/r','sin(2L\pi/C)');
  title(sprintf('L=%.0f',L));
  %hold off;
  shg;
  pause;
end
%%
Ra = (r^4 + [-1 1]*r^2*sqrt(r^4-rd^2*tand(th)^2*L^2))/(L*rd^2*tand(th)^2);
Rb = (4*r^2*L+[-1 1]*2*r*sqrt(4*r^2*L^2-L^2*W^2))/(W^2);
%%
L = 50:.2:100;
w = W/2;
r = W./(2*sin(2*L*pi/C));
% R = (r4 + r.^2.*sqrt(sq))./(L*rd^2*tand(th)^2);
R = (4*r.^2.*L+2*r.*L.*sqrt(4*r.^2-W^2))/(W^2);
plot(L,R); shg;
%%
plot(L,r); shg
%%
% This is the change in the axial coordinate between the center
% of a spherical optic and radius r.
dL = R-sqrt(R.^2-r.^2);
plot(L,dL); shg
%%
clf;
scales = 1:.1:1.5;
for scale = scales
  r = scale * W./(2*sin(2*L*pi/C));
  R = (4*r.^2.*L+2*r.*L.*sqrt(4*r.^2-W^2))/(W^2);
  plot(L,R);
  hold on;
end
hold off;
xlabel('Cavity Length');
ylabel('Radius of Curvature');
legend(num2str(scales')); shg;

%%
% Based on online chat with Edmund sales, radius of curvature is probably
% constrained by the ability to accurately measure the sag, and they appear
% to be able to easily measure 0.1 mm, so going much below that seems
% inadvisable.
clf;
gs = [0.001 0.005 0.01];
for g = gs % minimum sag
  r = 1:0.1:(3*2.54/2);
  R = (g^2+r.^2)/(2*g);
  D = r*2/2.54;
  semilogy(D,R); hold on;
end
hold off;
legend(num2str(gs'),'location','northwest');
xlabel('Mirror Diameter, inches');
ylabel('Max Radius of Curvature, cm');
title('RC vs Tolerance (cm)');
grid;
shg;

%%
% Another approach to constraints on R using overlap criteria
% Interleave==1,2,3a,3b
% %%% check w/r < rd tan(th) L/r^2
W = .4;
rd = .1;
th = 15;
C = 3000;
for L=50:10:110
  n1 = ceil(C/(2*L));
  p1 = pi/n1;
  w1 = W/2;
  r1 = w1/sin(p1);
  R1 = r1*L*(r1-[1,-1]*sqrt(r1^2-w1^2))/(w1^2);
  n2 = 1+2*ceil(C/(4*L)-1/2);
  p2 = 2*pi/n2;
  w2 = w1*sin(p2)/sin(pi/n2);
  r2 = w1/sin(pi/n2);
  R2 = r2*L*(r2-[1,-1]*sqrt(r2^2-w2^2))/(w2^2);
  n3a = ceil(C/(6*L)-1/3)*3+1;
  w3a = w1*sin(3*pi/n3a)/sin(pi/n3a);
  r3a = w1/sin(pi/n3a);
  R3a = r3a*L*(r3a-[1,-1]*sqrt(r3a^2-w3a^2))/(w3a^2);
  n3b = ceil(C/(6*L)-2/3)*3+2;
  w3b = w1*sin(3*pi/n3b)/sin(pi/n3b);
  r3b = w1/sin(pi/n3b);
  R3b = r3b*L*(r3b-[1,-1]*sqrt(r3b^2-w3b^2))/(w3b^2);
  p3a = 3*pi/n3a;
  p3b = 3*pi/n3b;
  R = (L/2):.1:500;
  if R1(end) > R(end)
    R1(end) = R(end);
  end
  if R2(end) > R(end)
    R2(end) = R(end);
  end
  if R3a(end) > R(end)
    R3a(end) = R(end);
  end
  if R3b(end) > R(end)
    R3b(end) = R(end);
  end
  
  w_r = (L./R).*sqrt(2*R/L - 1);
  Ra = L./(1+[1,-1]*cos(p1));
  clf;
  x = [R(1) R(end)];
  y = [1 1];
  plot(R,w_r,R1,sin(p1)*y,R2,sin(p2)*y,R3a,sin(p3a)*y,R3b,sin(p3b)*y);
  xlabel('R, cm');
  ylabel('w/r');
  title(sprintf('L=%.0f  n1=%d n2=%d n3a=%d n3b=%d', L, n1, n2, n3a, n3b));
  grid;
  legend('w/r','I=1','I=2','I=3a','I=3b');
  shg;
  % set(gca,'xlim',[L/2 L/2+5]);
  pause;
end
%%
L=50;
n3a = ceil(C/(6*L)-1/3)*3+1;
n3b = ceil(C/(6*L)-2/3)*3+2;
n3 = min(n3a,n3b);
p3 = 6*pi/n3;

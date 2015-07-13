%%
P.R1 = 75;
P.r1 = 3;
P.L = 50;
P.r2 = [0.5:.01:1.8];
P.Rw1 = 0.2;
R = exparam(P);
%% Radius of curvature for plano-convex 1" optics
ZCPXR = -0.1 * [
  71.46
  89.14
  106.15
  139.96
  180.32
  209.91
  281.20
  356.51
  698.17
];
ZCPXC = [1./ZCPXR; 0];
%%
clf;
np = 5;
ax = zeros(np,1);
for i=1:np
  ax(i) = nsubplot(np,1,i);
  switch i
    case 1
      plot(ax(i), R.r2, R.C2);
      ylabel(ax(i), 'IM2 Curvature, cm^{-1}');
      set(ax(i),'ytick',ZCPXC,'ygrid','on');
    case 2
      plot(ax(i), R.r2, R.r15);
      ylabel(ax(i),'radius at 15^o, cm');
      set(ax(i),'ytick',[0.1 0.2 0.5 1],'ygrid','on');
    case 3
      plot(ax(i), R.r2, R.RL);
      ylabel(ax(i),'RIM Space cm');
    case 4
      plot(ax(i), R.r2,R.Rr1*2/2.54);
      ylabel(ax(i),'RIM minimum diameter inches');
      set(ax(i),'ytick',[1 2 3 6 8],'ygrid','on');
    case 5
      plot(ax(i), R.r2,R.RR1/(2*2.54));
      ylabel(ax(i),'RIM Focal Length, inches');
      set(ax(i),'ytick',[4 8 12 24],'ygrid','on');
  end
end
title(ax(1),sprintf('IR1=%.0f IL=%.0f Ir1=%.1f', R.R1, R.L, R.r1));
set(ax(1:end-1),'xticklabel',[]);
set(ax(2:2:end),'yaxislocation','right');
xlabel(ax(np),'Ir2, cm');
linkaxes(ax,'x');
shg;
%%
figure;
P2 = P;
P2.r1 = 2;
R2 = exparam(P2);
plot(R.r15,R.RL,R2.r15,R2.RL);
legend('r1=3','r1=2');
xlabel('Radius at 15^o, cm');
ylabel('RIM Length, cm');
title(sprintf('IR1=%.0f IL=%.0f', R.R1, R.L));

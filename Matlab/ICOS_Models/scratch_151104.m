%%
% Investigate different interleave patterns and evaluate for relevant
% design criteria:
%   r_min = r1_min (min pattern radius to avoid overlap)
%   r_max = r1_max (max pattern radius to achieve focus)
C = 3000;
B = .4;
rd = .1;
th = 15;
L = 35;
m_min = ceil(C/(2*L)); % This spot can overlap
rdtanth = rd*tand(th);
k_set = 1:floor(m_min/2);
Summary.RLmin = zeros(size(k_set));
Summary.r_max = zeros(size(k_set));
Summary.R1 = zeros(size(k_set));
Summary.R2 = zeros(size(k_set));
Summary.m = zeros(size(k_set));
for i = 1:length(k_set);
  k = k_set(i);
  m = m_min;
  while gcd(m,k) ~= 1
    m = m+1;
  end
  phi = pi/m;
  Phi = k*phi;
  R1pole = L/sin(Phi)^2;
  if R1pole > 2000
    R1 = interp1([0 2000],[L 2000],1:2000)';
    R2 = L*(R1-L)./(sin(Phi)^2*R1-L);
    RR = []; % L/(1+cos(Phi));
    Rmax = 2000;
  else
    R1a = interp1([0 1001],[L R1pole],1:1000)';
    R1b = interp1([0 1000],[R1pole 2000],1:1000)';
    R1 = [R1a;R1pole;R1b];
    R2 = [L*(R1a-L)./(sin(Phi)^2*R1a-L); NaN; L*(R1b-L)./(sin(Phi)^2*R1b-L)];
    RR = L./(1+cos(Phi)*[-1]);
    Rmax = min(RR*1.1,2000);
    RR = RR(RR<2000);
  end
  R2(R2 > 2000 | R2 < -2000) = NaN;
  % fprintf(1,'Break at R1 = %.1f\n', L/sin(Phi)^2);
  r_min = B/(2*sin(pi/m));
  r_max = sqrt(rdtanth*sqrt(((R2-L).*R1.^2*L)./((R1-L).*(R1+R2-L))));
  rRR = interp1(R1,r_max,RR,'linear','extrap');
% This (commented) calculation is wrong, but I haven't figured out why yet
%   k2 = ((R1-L).*(R1.*R2-L))./((R2-L).*R1.^2.*L);
%   k2(k2<0) = NaN;
%   s1 = r_max.*sqrt(k2);
  r2 = abs(r_max.*cos(Phi).*R2./(R2-L));
  w2 = abs(r_max*sin(Phi)*cos(Phi).*R2./(R2-L));
  s1 = w2/L;

  clf;
  ax = [nsubplot(3,1,1) nsubplot(3,1,2) nsubplot(3,1,3)];
  plot(ax(1),R1,r_max,'b',[R1(1) R1(end)],r_min*[1 1],'g',RR,rRR,'*r',...
    R1,r2);
  % xlabel(ax(1),'R1 cm');
  ylabel(ax(1),'r cm');
  xlim(ax(1),[L/2 Rmax]);
  ylim(ax(1),[0 5]);
  set(ax(1),'XTickLabel',[],'YAxisLocation','Right');
  title(ax(1),sprintf('L=%d %d/%d', L, k, m));
  
  plot(ax(3),R1,R2);
  ylabel(ax(3),'R_2');
  set(ax(3),'YAxisLocation','Right');
  
  Rw1 = 1.1*B/2;
  RL = abs(Rw1./s1);
  rRL = interp1(R1,RL,RR,'linear','extrap');
  vok = r_max > r_min & r_max >= r2;
  vnok = ~vok;
  RLok = RL;
  RLok(vnok) = NaN;
  RLnok = RL;
  RLnok(vok) = NaN;
  
  Summary.RLmin(i) = nanmin(RLok);
  mi = find(RLok == Summary.RLmin(i));
  Summary.R1(i) = R1(mi);
  Summary.R2(i) = R2(mi);
  Summary.r_max(i) = r_max(mi);
  Summary.m(i) = m;
  
  plot(ax(2),R1,RLok,'b',R1,RLnok,'r',RR,rRL,'*r');
  xlim(ax(2),[L/2 Rmax]);
  ylim(ax(2),[0 min(L,nanmax(RL))]);
  xlabel(ax(end),'R1 cm');
  ylabel(ax(2),'RL cm');
  set(ax(2),'XTickLabel',[],'YAxisLocation','Left');
  % title(sprintf('RL: L=%d %d/%d', L, k, m));
  linkaxes(ax,'x');
  
  ii = find(~isnan(RLok),1,'last');
  SP.R1 = R1(ii);
  SP.R2 = R2(ii);
  SP.L = L;
  SP.Rw1 = 1.1*B/2;
  SP.RL = SP.Rw1/s1(ii);
  Res = exparam(SP);
  check_params(k, Res);
  pause;
end
%%
clf;
ax = [ nsubplot(3,1,1), nsubplot(3,1,2), nsubplot(3,1,3)];
plot(ax(1),k_set, Summary.RLmin, '*');
plot(ax(2),k_set, Summary.r_max, '*');
plot(ax(3),k_set, Summary.R1, '*', k_set, Summary.R2, '+');
set(ax(1),'XTickLabel',[]);
set(ax(2),'XTickLabel',[],'YAxisLocation','Right');
ylabel(ax(1),'RL_{min} cm');
ylabel(ax(2),'r_{1,max} cm');
ylabel(ax(3),'R_{1,2}');
xlabel(ax(3),'k');
linkaxes(ax, 'x');
xlim(ax(3),minmax(k_set));
%%
figure;
% v = Summary.R1 < 200 & Summary.R2 > -200;
v = Summary.R1 > 0;
scatter(Summary.R1(v),Summary.R2(v),[],Summary.m(v));
xl = xlim;
R1max = xl(2);
hold on;
plot([0, R1max], [0, -R1max],'k');
colorbar;


%%
% Demonstrate different interleave patterns
C = 3000;
B = .4;
rd = .1;
th = 15;
L = 35;
m_min = ceil(C/(2*L)); % This spot can overlap
rdtanth = rd*tand(th);
for k=1:floor(m_min/2)
  m = m_min;
  while gcd(m,k) ~= 1
    m = m+1;
  end
  phi = pi/m;
  Phi = k*phi;
  R1pole = L/sin(Phi)^2;
  if R1pole > 2000
    R1 = interp1([0 2000],[L/2 2000],1:2000)';
    R2 = L*(R1-L)./(sin(Phi)^2*R1-L);
    RR = L/(1+cos(Phi));
    Rmax = 2000;
  else
    R1a = interp1([0 1001],[L/2 R1pole],1:1000)';
    R1b = interp1([0 1000],[R1pole 2000],1:1000)';
    R1 = [R1a;R1pole;R1b];
    R2 = [L*(R1a-L)./(sin(Phi)^2*R1a-L); NaN; L*(R1b-L)./(sin(Phi)^2*R1b-L)];
    RR = L./(1+cos(Phi)*[1,-1]);
    Rmax = min(RR(2)*1.1,2000);
    RR = RR(RR<2000);
  end
  R2(R2 > 2000) = NaN;
  % fprintf(1,'Break at R1 = %.1f\n', L/sin(Phi)^2);
  r_min = B/(2*sin(pi/m));
  r_max = sqrt(rdtanth*sqrt(((R2-L).*R1.^2*L)./((R1-L).*(R1+R2-L))));
  rRR = interp1(R1,r_max,RR,'linear','extrap');
  clf;
  ax = [nsubplot(3,1,1) nsubplot(3,1,2) nsubplot(3,1,3)];
  plot(ax(1),R1,r_max,'b',[R1(1) R1(end)],r_min*[1 1],'g',RR,rRR,'*r');
  % xlabel(ax(1),'R1 cm');
  ylabel(ax(1),'r_{max} cm');
  xlim(ax(1),[L/2 Rmax]);
  ylim(ax(1),[0 5]);
  set(ax(1),'XTickLabel',[],'YAxisLocation','Right');
  title(ax(1),sprintf('L=%d %d/%d', L, k, m));
  
%   k2 = ((R1-L).*(R1.*R2-L))./((R2-L).*R1.^2.*L);
%   k2(k2<0) = NaN;
%   s1 = r_max.*sqrt(k2);
  w2 = r_max*sin(Phi)*cos(Phi).*R2./(R2-L);
  s1 = w2/L;
  plot(ax(3),R1,s1);
  ylabel(ax(3),'s1');
  set(ax(3),'YAxisLocation','Right');
  
  Rw1 = 1.1*B/2;
  RL = abs(Rw1./s1);
  rRL = interp1(R1,RL,RR,'linear','extrap');
  vok = r_max > r_min;
  vnok = ~vok;
  RLok = RL;
  RLok(vnok) = NaN;
  RLnok = RL;
  RLnok(vok) = NaN;
  plot(ax(2),R1,RLok,'b',R1,RLnok,'r',RR,rRL,'*r');
  xlim(ax(2),[L/2 Rmax]);
  ylim(ax(2),[0 min(L,nanmax(RL))]);
  xlabel(ax(end),'R1 cm');
  ylabel(ax(2),'RL cm');
  set(ax(2),'XTickLabel',[],'YAxisLocation','Left');
  % title(sprintf('RL: L=%d %d/%d', L, k, m));
  linkaxes(ax,'x');
  pause;
end
